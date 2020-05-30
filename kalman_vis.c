/*
 * CDDL HEADER START
 *
 * The contents of this file are subject to the terms of the
 * Common Development and Distribution License, Version 1.0 only
 * (the "License").  You may not use this file except in compliance
 * with the License.
 *
 * You can obtain a copy of the license in the file COPYING
 * or http://www.opensource.org/licenses/CDDL-1.0.
 * See the License for the specific language governing permissions
 * and limitations under the License.
 *
 * When distributing Covered Code, include this CDDL HEADER in each
 * file and include the License file COPYING.
 * If applicable, add the following below this CDDL HEADER, with the
 * fields enclosed by brackets "[]" replaced with your own identifying
 * information: Portions Copyright [yyyy] [name of copyright owner]
 *
 * CDDL HEADER END
 */
/*
 * Copyright 2020 Saso Kiselkov. All rights reserved.
 */

#include <XPLMDisplay.h>

#include <acfutils/geom.h>
#include <acfutils/list.h>
#include <acfutils/math.h>
#include <acfutils/mt_cairo_render.h>
#include <acfutils/safe_alloc.h>
#include <acfutils/widget.h>
#include <acfutils/thread.h>

#include "kalman_assert.h"
#include "kalman_vis.h"

#define	GRAPH_WIDTH		600
#define	GRAPH_HEIGHT		120
#define	GRAPH_DATA_WIDTH	80
#define	COV_COLUMN_WIDTH	120
#define	RENDER_FPS		20
#define	MAX_SAMPLES		200
#define	PX_PER_SAMPLE		3

typedef struct {
	kalman_vec_t		m;
	kalman_vec_t		state;
	list_node_t		node;
} sample_t;

struct kalman_vis_s {
	kalman_t		*kal;		/* immutable */
	unsigned		state_len;	/* immutable */
	mt_cairo_render_t	*mtcr;		/* immutable */
	XPLMWindowID		win;		/* immutable */
	win_resize_ctl_t	winctl;		/* immutable */

	mutex_t			lock;
	/* protected by lock */
	list_t			samples;
	kalman_mat_t		cov;
	char			labels[KALMAN_VEC_LEN][128];

	unsigned		decimals[KALMAN_VEC_LEN];	/* atomic */
	unsigned		cov_precision;	/* atomic */
};

static const vect3_t color_table[] = {
    {1, 0, 0},
    {0, 0.8, 0},
    {0, 0, 1},
    {1, 0, 0},
    {0, 0.8, 0},
    {0, 0, 1},
    {1, 1, 0},
    {1, 0, 1},
    {0, 1, 1},
    {0.8, 0.5, 0},
    {0.8, 0, 0.5}
};

static void
win_draw(XPLMWindowID win, void *refcon)
{
	kalman_vis_t *vis;
	int left, top, right, bottom, w, h;

	KAL_ASSERT(win != NULL);
	KAL_ASSERT(refcon != NULL);
	vis = refcon;
	win_resize_ctl_update(&vis->winctl);
	XPLMGetWindowGeometry(win, &left, &top, &right, &bottom);
	w = right - left;
	h = top - bottom;

	KAL_ASSERT(vis->mtcr != NULL);
	mt_cairo_render_draw(vis->mtcr, VECT2(left, bottom), VECT2(w, h));
}

static void
render_graph(cairo_t *cr, kalman_vis_t *vis, unsigned idx)
{
	char buf[32];
	int i;
	vect3_t color;
	kalman_real_t minval = INFINITY, maxval = -INFINITY;
	kalman_real_t stateval;
	cairo_text_extents_t te;

	KAL_ASSERT(cr != NULL);
	KAL_ASSERT(vis != NULL);
	KAL_ASSERT3U(idx, <, vis->state_len);
	KAL_ASSERT(list_count(&vis->samples) != 0);

	cairo_save(cr);
	cairo_translate(cr, GRAPH_DATA_WIDTH, (idx + 0.5) * GRAPH_HEIGHT);
	cairo_rectangle(cr, 0.5, -GRAPH_HEIGHT / 2 + 0.5,
	    GRAPH_WIDTH - 1, GRAPH_HEIGHT - 1);
	cairo_rectangle(cr, -GRAPH_DATA_WIDTH + 0.5, -GRAPH_HEIGHT / 2 + 0.5,
	    GRAPH_DATA_WIDTH - 1, GRAPH_HEIGHT - 1);
	cairo_stroke(cr);

	/* Determine min/max values for this graph */
	for (sample_t *sample = list_head(&vis->samples); sample != NULL;
	    sample = list_next(&vis->samples, sample)) {
		minval = MIN(minval, sample->m.v[idx]);
		minval = MIN(minval, sample->state.v[idx]);
		maxval = MAX(maxval, sample->m.v[idx]);
		maxval = MAX(maxval, sample->state.v[idx]);
	}
	stateval = ((sample_t *)list_head(&vis->samples))->state.v[idx];
	ASSERT(isfinite(minval));
	ASSERT(isfinite(maxval));

	/* Draw the measurement values - those are all black */
	i = 0;
	for (sample_t *sample = list_head(&vis->samples); sample != NULL;
	    sample = list_next(&vis->samples, sample), i++) {
		kalman_real_t val = sample->m.v[idx];
		double x = GRAPH_WIDTH - i * PX_PER_SAMPLE;
		double y = round(fx_lin(val, minval, GRAPH_HEIGHT / 2,
		    maxval, -GRAPH_HEIGHT / 2));

		cairo_rectangle(cr, x - 1, y - 1, 2, 2);
	}
	cairo_fill(cr);

	/* Draw filter state */
	color = color_table[idx % ARRAY_NUM_ELEM(color_table)];
	cairo_set_source_rgb(cr, color.x, color.y, color.z);
	i = 0;
	for (sample_t *sample = list_head(&vis->samples); sample != NULL;
	    sample = list_next(&vis->samples, sample), i++) {
		kalman_real_t val = sample->state.v[idx];
		double x = GRAPH_WIDTH - i * PX_PER_SAMPLE;
		double y = fx_lin(val, minval, GRAPH_HEIGHT / 2,
		    maxval, -GRAPH_HEIGHT / 2);

		if (i == 0)
			cairo_move_to(cr, x + 0.5, y + 0.5);
		else
			cairo_line_to(cr, x + 0.5, y + 0.5);
	}
	cairo_stroke(cr);

	/* Draw the state's current value */
	snprintf(buf, sizeof (buf), "%.*f",
	    fixed_decimals(stateval, vis->decimals[idx]),
	    (double)stateval);
	cairo_text_extents(cr, buf, &te);
	cairo_move_to(cr, -GRAPH_DATA_WIDTH / 2 - te.width / 2,
	    -te.height / 2 - te.y_bearing);
	cairo_show_text(cr, buf);

	/* Draw the minimum and maximum graph values */
	cairo_set_source_rgb(cr, 0, 0, 0);
	snprintf(buf, sizeof (buf), "%.*f",
	    fixed_decimals(maxval, vis->decimals[idx]), (double)maxval);
	cairo_text_extents(cr, buf, &te);
	cairo_move_to(cr, -GRAPH_DATA_WIDTH / 2 - te.width / 2,
	    -GRAPH_HEIGHT / 2 + te.height - te.y_bearing);
	cairo_show_text(cr, buf);

	snprintf(buf, sizeof (buf), "%.*f",
	    fixed_decimals(minval, vis->decimals[idx]), (double)minval);
	cairo_text_extents(cr, buf, &te);
	cairo_move_to(cr, -GRAPH_DATA_WIDTH / 2 - te.width / 2,
	    GRAPH_HEIGHT / 2 - te.height - te.y_bearing - 10);
	cairo_show_text(cr, buf);

	/* Draw the label */
	if (strlen(vis->labels[idx]) != 0) {
		enum { MARGIN = 2 };

		cairo_text_extents(cr, vis->labels[idx], &te);

		cairo_set_source_rgb(cr, 0.85, 0.85, 0.85);
		cairo_rectangle(cr, 1, -GRAPH_HEIGHT / 2 + 1,
		    te.width + 2 * MARGIN, te.height + 2 * MARGIN);
		cairo_fill(cr);

		cairo_set_source_rgb(cr, 0, 0, 0);
		cairo_rectangle(cr, 0.5, -GRAPH_HEIGHT / 2 + 0.5,
		    te.width + 2 * MARGIN, te.height + 2 * MARGIN);
		cairo_stroke(cr);

		cairo_set_source_rgb(cr, 0, 0, 0);
		cairo_move_to(cr, MARGIN, -GRAPH_HEIGHT / 2 + MARGIN -
		    te.y_bearing);
		cairo_show_text(cr, vis->labels[idx]);
	}

	cairo_restore(cr);
}

static void
render_cov(cairo_t *cr, kalman_vis_t *vis)
{
	cairo_text_extents_t te;

	KAL_ASSERT(cr != NULL);
	KAL_ASSERT(vis != NULL);

	cairo_save(cr);
	cairo_translate(cr, GRAPH_DATA_WIDTH + GRAPH_WIDTH, 0);

	cairo_text_extents(cr, "Covariance", &te);
	cairo_move_to(cr,
	    (vis->state_len * COV_COLUMN_WIDTH) / 2 - te.width / 2, 20);
	cairo_show_text(cr, "Covariance");

	for (unsigned x = 0; x < vis->state_len; x++) {
		for (unsigned y = 0; y < vis->state_len; y++) {
			char buf[32];

			snprintf(buf, sizeof (buf), "%.*f", vis->cov_precision,
			    (double)KALMAN_MATxy(vis->cov, x, y));
			cairo_text_extents(cr, buf, &te);
			cairo_move_to(cr,
			    (x + 0.5) * COV_COLUMN_WIDTH - te.width / 2,
			    (y + 0.5) * COV_COLUMN_WIDTH - te.height / 2 -
			    te.y_bearing);
			cairo_show_text(cr, buf);
		}
	}

	cairo_restore(cr);
}

static void
kal_vis_render(cairo_t *cr, unsigned w, unsigned h, void *userinfo)
{
	kalman_vis_t *vis;

	KAL_ASSERT(cr != NULL);
	KAL_ASSERT(userinfo != NULL);
	vis = userinfo;

	cairo_set_source_rgb(cr, 1, 1, 1);
	cairo_paint(cr);

	cairo_set_source_rgb(cr, 0, 0, 0);
	cairo_set_font_size(cr, 20);
	cairo_set_line_width(cr, 1);

	mutex_enter(&vis->lock);

	if (list_count(&vis->samples) != 0) {
		for (unsigned i = 0; i < vis->state_len; i++)
			render_graph(cr, vis, i);
		render_cov(cr, vis);
	} else {
		cairo_text_extents_t te;
		cairo_set_source_rgb(cr, 1, 0, 0);
		cairo_text_extents(cr, "NO DATA", &te);
		cairo_move_to(cr, w / 2 - te.width / 2,
		    h / 2 - te.height / 2 - te.y_bearing);
		cairo_show_text(cr, "NO DATA");
	}

	mutex_exit(&vis->lock);
}

kalman_vis_t *
kalman_vis_alloc(kalman_t *kal, const char *name)
{
	kalman_vis_t *vis = safe_calloc(1, sizeof (*vis));
	XPLMCreateWindow_t cr = {
	    .structSize = sizeof (cr),
	    .left = 100,
	    .bottom = 100,
	    .drawWindowFunc = win_draw,
	    .refcon = vis,
	    .decorateAsFloatingWindow = xplm_WindowDecorationRoundRectangle,
	    .layer = xplm_WindowLayerFloatingWindows
	};
	unsigned w, h;

	KAL_ASSERT(kal != NULL);
	vis->state_len = kalman_get_state_len(kal);
	KAL_ASSERT3U(vis->state_len, <, KALMAN_VEC_LEN);
	KAL_ASSERT(name != NULL);
	for (unsigned i = 0; i < vis->state_len; i++)
		vis->decimals[i] = 8;
	vis->cov_precision = 3;

	w = GRAPH_WIDTH + GRAPH_DATA_WIDTH + vis->state_len * COV_COLUMN_WIDTH;
	h = vis->state_len * GRAPH_HEIGHT;
	cr.top = 100 + h;
	cr.right = 100 + w;

	vis->kal = kal;
	mutex_init(&vis->lock);
	list_create(&vis->samples, sizeof (sample_t), offsetof(sample_t, node));
	vis->mtcr = mt_cairo_render_init(w, h, RENDER_FPS, NULL,
	    kal_vis_render, NULL, vis);
	vis->win = XPLMCreateWindowEx(&cr);
	KAL_ASSERT(vis->win != NULL);
	classic_win_center(vis->win);
	win_resize_ctl_init(&vis->winctl, vis->win, w, h);
	XPLMSetWindowTitle(vis->win, name);

	return (vis);
}

void
kalman_vis_free(kalman_vis_t *vis)
{
	sample_t *sample;

	if (vis == NULL)
		return;

	mt_cairo_render_fini(vis->mtcr);
	while ((sample = list_remove_head(&vis->samples)) != NULL)
		free(sample);
	list_destroy(&vis->samples);
	mutex_destroy(&vis->lock);
	XPLMDestroyWindow(vis->win);
	memset(vis, 0, sizeof (*vis));
	free(vis);
}

void
kalman_vis_update(kalman_vis_t *vis, const kalman_vec_t *m)
{
	sample_t *sample = safe_calloc(1, sizeof (*sample));

	KAL_ASSERT(vis != NULL);
	KAL_ASSERT(m != NULL);

	sample->m = *m;
	sample->state = kalman_get_state(vis->kal);

	mutex_enter(&vis->lock);
	while (list_count(&vis->samples) >= MAX_SAMPLES) {
		sample_t *old_sample = list_remove_tail(&vis->samples);
		free(old_sample);
	}
	list_insert_head(&vis->samples, sample);
	vis->cov = kalman_get_cov_mat(vis->kal);
	mutex_exit(&vis->lock);
}

void
kalman_vis_reset(kalman_vis_t *vis)
{
	sample_t *sample;

	KAL_ASSERT(vis != NULL);

	mutex_enter(&vis->lock);
	while ((sample = list_remove_head(&vis->samples)) != NULL)
		free(sample);
	mutex_exit(&vis->lock);
}

void
kalman_vis_open(kalman_vis_t *vis)
{
	KAL_ASSERT(vis != NULL);
	KAL_ASSERT(vis->win != NULL);
	XPLMSetWindowIsVisible(vis->win, true);
	XPLMBringWindowToFront(vis->win);
}

bool
kalman_vis_is_open(const kalman_vis_t *vis)
{
	KAL_ASSERT(vis != NULL);
	KAL_ASSERT(vis->win != NULL);
	return (XPLMGetWindowIsVisible(vis->win));
}

void
kalman_vis_set_decimals(kalman_vis_t *vis, unsigned state_var,
    unsigned decimals)
{
	KAL_ASSERT(vis != NULL);
	KAL_ASSERT3U(state_var, <, vis->state_len);
	vis->decimals[state_var] = decimals;
}

void
kalman_vis_set_label(kalman_vis_t *vis, unsigned state_var, const char *label)
{
	KAL_ASSERT(vis != NULL);
	KAL_ASSERT3U(state_var, <, vis->state_len);
	KAL_ASSERT(label != NULL);
	mutex_enter(&vis->lock);
	strlcpy(vis->labels[state_var], label, sizeof (vis->labels[state_var]));
	mutex_exit(&vis->lock);
}

void
kalman_vis_set_cov_precision(kalman_vis_t *vis, unsigned decimals)
{
	KAL_ASSERT(vis != NULL);
	vis->cov_precision = decimals;
}
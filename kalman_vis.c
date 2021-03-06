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

#include <stdarg.h>

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

#define	GRAPH_WIDTH		(vis->max_samples * vis->px_per_sample)
#define	GRAPH_DATA_WIDTH	110
#define	COV_COLUMN		100
#define	CONT_COLUMN		100
#define	RENDER_FPS		20

#define	HEADING_FONT_SZ		21
#define	COV_DATA_FONT_SZ	16
#define	GRAPH_DATA_FONT_SZ	18

#define	AUTO_DECIMALS(value, decimals) \
	((decimals) < 0 ? -(decimals) : fixed_decimals((value), (decimals)))

typedef struct {
	kalman_dmat_t		*m;
	kalman_dmat_t		*state;
	list_node_t		node;
} sample_t;

struct kalman_vis_s {
	/* immutable */
	kalman_t		*kal;
	unsigned		state_len;
	unsigned		max_samples;
	double			px_per_sample;
	double			row_height;
	mt_cairo_render_t	*mtcr;
	XPLMWindowID		win;
	win_resize_ctl_t	winctl;

	mutex_t			lock;
	/* protected by lock */
	list_t			samples;
	kalman_dmat_t		*cov;
	kalman_dmat_t		*m_cov;
	uint64_t		m_cov_age;
	kalman_dmat_t		*cont;
	char			labels[KALMAN_VEC_LEN][128];

	int			decimals[KALMAN_VEC_LEN];	/* atomic */
	int			cov_precision;	/* atomic */
};

static const vect3_t color_table[] = {
	{1,	0,	0},	/* Red */
	{0,	0.67,	0},	/* Green */
	{0,	0,	1},	/* Blue */
	{0.67,	0,	0.67},	/* Magenta */
	{0.67,	0.33,	0},	/* Brown */
	{0,	0.67,	0.67},	/* Cyan */
	{0.33,	0,	0.67},	/* Purple */
	{0.9,	0,	0.47},	/* Pink */
	{0.67,	0.67,	0},	/* Yellow */
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
render_centered_text(cairo_t *cr, double x, double y, const char *format, ...)
{
	char buf[128];
	va_list ap;
	cairo_text_extents_t te;

	va_start(ap, format);
	vsnprintf(buf, sizeof (buf), format, ap);
	va_end(ap);

	cairo_text_extents(cr, buf, &te);
	cairo_move_to(cr, x - te.width / 2, y - te.height / 2 - te.y_bearing);
	cairo_show_text(cr, buf);
}

static void
render_graph(cairo_t *cr, kalman_vis_t *vis, unsigned idx)
{
	int i;
	vect3_t color;
	kalman_real_t minval = INFINITY, maxval = -INFINITY;
	kalman_real_t stateval;
	sample_t *sample;

	KAL_ASSERT(cr != NULL);
	KAL_ASSERT(vis != NULL);
	KAL_ASSERT3U(idx, <, vis->state_len);
	KAL_ASSERT(list_count(&vis->samples) != 0);

	cairo_save(cr);
	cairo_translate(cr, GRAPH_DATA_WIDTH, (idx + 0.5) * vis->row_height);
	cairo_rectangle(cr, 0.5, -vis->row_height / 2 + 0.5,
	    GRAPH_WIDTH - 1, vis->row_height - 1);
	cairo_rectangle(cr, -GRAPH_DATA_WIDTH + 0.5, -vis->row_height / 2 + 0.5,
	    GRAPH_DATA_WIDTH - 1, vis->row_height - 1);
	cairo_stroke(cr);

	/* Determine min/max values for this graph */
	for (sample = list_head(&vis->samples); sample != NULL;
	    sample = list_next(&vis->samples, sample)) {
		if (sample->m != NULL &&
		    !isnan(KAL_DMATyx(*sample->m, idx, 0))) {
			minval = MIN(minval, KAL_DMATyx(*sample->m, idx, 0));
			maxval = MAX(maxval, KAL_DMATyx(*sample->m, idx, 0));
		}
		ASSERT(sample->state != NULL);
		minval = MIN(minval, KAL_DMATyx(*sample->state, idx, 0));
		maxval = MAX(maxval, KAL_DMATyx(*sample->state, idx, 0));
	}
	sample = list_head(&vis->samples);
	ASSERT(sample != NULL);
	stateval = KAL_DMATyx(*sample->state, idx, 0);
	ASSERT(isfinite(minval));
	ASSERT(isfinite(maxval));

	/* Draw the measurement values - those are all black */
	for (sample = list_head(&vis->samples), i = 0; sample != NULL;
	    sample = list_next(&vis->samples, sample), i++) {
		double val, x, y;

		if (sample->m == NULL || isnan(KAL_DMATyx(*sample->m, idx, 0)))
			continue;
		val = KAL_DMATyx(*sample->m, idx, 0);
		x = GRAPH_WIDTH - i * vis->px_per_sample;
		y = fx_lin(val, minval, vis->row_height / 2,
		    MAX(maxval, minval + 1e-6), -vis->row_height / 2);

		cairo_rectangle(cr, x - 1, y - 1, 2, 2);
	}
	cairo_fill(cr);

	/* Draw filter state */
	color = color_table[idx % ARRAY_NUM_ELEM(color_table)];
	cairo_set_source_rgb(cr, color.x, color.y, color.z);
	for (sample = list_head(&vis->samples), i = 0; sample != NULL;
	    sample = list_next(&vis->samples, sample), i++) {
		double val = KAL_DMATyx(*sample->state, idx, 0);
		double x = GRAPH_WIDTH - i * vis->px_per_sample;
		double y = fx_lin(val, minval, vis->row_height / 2,
		    MAX(maxval, minval + 1e-6), -vis->row_height / 2);

		ASSERT(!isnan(val));
		if (i == 0)
			cairo_move_to(cr, x + 0.5, y + 0.5);
		else
			cairo_line_to(cr, x + 0.5, y + 0.5);
	}
	cairo_stroke(cr);

	/* Draw the state's current value */
	render_centered_text(cr, -GRAPH_DATA_WIDTH / 2, 0, "%.*f",
	    AUTO_DECIMALS(stateval, vis->decimals[idx]), (double)stateval);

	/* Draw the minimum and maximum graph values */
	cairo_set_source_rgb(cr, 0, 0, 0);
	render_centered_text(cr, -GRAPH_DATA_WIDTH / 2,
	    -vis->row_height / 2 + 15, "%.*f",
	    AUTO_DECIMALS(maxval, vis->decimals[idx]), (double)maxval);
	render_centered_text(cr, -GRAPH_DATA_WIDTH / 2,
	    vis->row_height / 2 - 15, "%.*f",
	    AUTO_DECIMALS(minval, vis->decimals[idx]), (double)minval);

	/* Draw the label */
	if (strlen(vis->labels[idx]) != 0) {
		enum { MARGIN = 2 };
		cairo_text_extents_t te;

		cairo_text_extents(cr, vis->labels[idx], &te);

		cairo_set_source_rgb(cr, 0.85, 0.85, 0.85);
		cairo_rectangle(cr, 1, -vis->row_height / 2 + 1,
		    te.width + 2 * MARGIN, te.height + 2 * MARGIN);
		cairo_fill(cr);

		cairo_set_source_rgb(cr, 0, 0, 0);
		cairo_rectangle(cr, 0.5, -vis->row_height / 2 + 0.5,
		    te.width + 2 * MARGIN, te.height + 2 * MARGIN);
		cairo_stroke(cr);

		cairo_set_source_rgb(cr, 0, 0, 0);
		cairo_move_to(cr, MARGIN, -vis->row_height / 2 + MARGIN -
		    te.y_bearing);
		cairo_show_text(cr, vis->labels[idx]);
	}

	cairo_restore(cr);
}

static void
render_cov(cairo_t *cr, kalman_vis_t *vis)
{
	KAL_ASSERT(cr != NULL);
	KAL_ASSERT(vis != NULL);
	KAL_ASSERT(vis->cov != NULL);

	cairo_save(cr);
	cairo_translate(cr, GRAPH_DATA_WIDTH + GRAPH_WIDTH, 0);

	cairo_set_font_size(cr, HEADING_FONT_SZ);
	render_centered_text(cr, (vis->state_len * COV_COLUMN) / 2, 10,
	    "Filter covariance");
	render_centered_text(cr, (vis->state_len * COV_COLUMN) / 2, 30,
	    (vis->state_len > 2 ? "(Measurement covariance)" :
	    "(Measurement cov.)"));

	cairo_set_font_size(cr, COV_DATA_FONT_SZ);
	for (unsigned x = 0; x < vis->cov->cols; x++) {
		for (unsigned y = 0; y < vis->cov->rows; y++) {
			double val;

			val = KAL_DMATyx(*vis->cov, y, x);
			if (val == 0) {
				render_centered_text(cr, (x + 0.5) * COV_COLUMN,
				    (y + 0.5) * vis->row_height, "0");
			} else {
				render_centered_text(cr, (x + 0.5) * COV_COLUMN,
				    (y + 0.5) * vis->row_height, "%.*f",
				    AUTO_DECIMALS(val, vis->cov_precision),
				    val);
			}

			if (vis->m_cov == NULL)
				continue;
			if (vis->m_cov_age > 1) {
				float grey = clampi(vis->m_cov_age, 0, 20) /
				    40.0;
				cairo_set_source_rgb(cr, grey, grey, grey);
			}
			val = KAL_DMATyx(*vis->m_cov, y, x);
			if (val == 0) {
				render_centered_text(cr, (x + 0.5) * COV_COLUMN,
				    (y + 0.5) * vis->row_height + 20, "(0)");
			} else if (!isnan(val)) {
				render_centered_text(cr, (x + 0.5) * COV_COLUMN,
				    (y + 0.5) * vis->row_height + 20, "(%.*f)",
				    AUTO_DECIMALS(val, vis->cov_precision),
				    val);
			} else {
				render_centered_text(cr, (x + 0.5) * COV_COLUMN,
				    (y + 0.5) * vis->row_height + 20, "(null)");
			}
			if (vis->m_cov_age > 1)
				cairo_set_source_rgb(cr, 0, 0, 0);
		}
	}

	cairo_restore(cr);
}

static void
render_cont(cairo_t *cr, unsigned h, kalman_vis_t *vis)
{
	KAL_ASSERT(cr != NULL);
	KAL_ASSERT(vis != NULL);
	KAL_ASSERT(vis->cont != NULL);

	cairo_save(cr);
	cairo_translate(cr, GRAPH_DATA_WIDTH + GRAPH_WIDTH +
	    vis->state_len * COV_COLUMN, 0);

	cairo_set_line_width(cr, 1);
	cairo_move_to(cr, 0.5, 0);
	cairo_line_to(cr, 0.5, h);
	cairo_stroke(cr);

	cairo_set_font_size(cr, HEADING_FONT_SZ);
	render_centered_text(cr, COV_COLUMN / 2, 10, "Control");
	render_centered_text(cr, COV_COLUMN / 2, 30, "vector");

	cairo_set_font_size(cr, COV_DATA_FONT_SZ);
	for (unsigned i = 0; i < vis->state_len; i++) {
		double val = KAL_DMATyx(*vis->cont, i, 0);

		if (val == 0) {
			render_centered_text(cr, CONT_COLUMN / 2,
			    (i + 0.5) * vis->row_height, "0");
		} else {
			render_centered_text(cr, CONT_COLUMN / 2,
			    (i + 0.5) * vis->row_height, "%.*f",
			    AUTO_DECIMALS(val, vis->cov_precision), val);
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
	cairo_set_font_size(cr, GRAPH_DATA_FONT_SZ);
	cairo_set_line_width(cr, 1);

	mutex_enter(&vis->lock);

	if (list_count(&vis->samples) != 0) {
		for (unsigned i = 0; i < vis->state_len; i++)
			render_graph(cr, vis, i);
		render_cov(cr, vis);
		render_cont(cr, h, vis);
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

static void
free_sample(sample_t *sample)
{
	KAL_ASSERT(sample != NULL);

	free(sample->m);
	free(sample->state);
	free(sample);
}

kalman_vis_t *
kalman_vis_alloc(kalman_t *kal, const char *name, unsigned max_samples,
    double px_per_sample, double row_height, mt_cairo_uploader_t *mtul)
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
	KAL_ASSERT(max_samples != 0);
	KAL_ASSERT(px_per_sample > 0);
	vis->state_len = kalman_get_state_len(kal);
	vis->max_samples = max_samples;
	vis->px_per_sample = px_per_sample;
	vis->row_height = row_height;
	KAL_ASSERT3U(vis->state_len, <=, KALMAN_VEC_LEN);
	KAL_ASSERT(name != NULL);
	for (unsigned i = 0; i < vis->state_len; i++)
		vis->decimals[i] = 8;
	vis->cov_precision = 7;

	w = GRAPH_WIDTH + GRAPH_DATA_WIDTH + vis->state_len * COV_COLUMN +
	    CONT_COLUMN;
	h = vis->state_len * vis->row_height;
	cr.top = 100 + h;
	cr.right = 100 + w;

	vis->kal = kal;
	vis->m_cov = NULL;
	mutex_init(&vis->lock);
	list_create(&vis->samples, sizeof (sample_t), offsetof(sample_t, node));
	vis->mtcr = mt_cairo_render_init(w, h, RENDER_FPS, NULL,
	    kal_vis_render, NULL, vis);
	ASSERT(vis->mtcr != NULL);
	mt_cairo_render_set_uploader(vis->mtcr, mtul);
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
		free_sample(sample);
	list_destroy(&vis->samples);
	mutex_destroy(&vis->lock);
	XPLMDestroyWindow(vis->win);
	memset(vis, 0, sizeof (*vis));
	free(vis);
}

void
kalman_vis_update(kalman_vis_t *vis, const kalman_vec_t *m,
    const kalman_mat_t *m_cov)
{
	kalman_dmat_t *dm, *dm_cov;

	if (m == NULL) {
		kalman_vis_dupdate(vis, NULL, NULL);
		return;
	}

	dm = kalman_dmat_alloc(vis->state_len, 1);
	dm_cov = kalman_dmat_alloc(vis->state_len, vis->state_len);
	for (unsigned row = 0; row < vis->state_len; row++) {
		KAL_DMATyx(*dm, row, 0) = KALMAN_VECi(*m, row);
		for (unsigned col = 0; col < vis->state_len; col++) {
			KAL_DMATyx(*dm_cov, row, col) =
			    KALMAN_MATyx(*m_cov, row, col);
		}
	}

	kalman_vis_dupdate(vis, dm, dm_cov);

	free(dm);
	free(dm_cov);
}

void
kalman_vis_dupdate(kalman_vis_t *vis, const kalman_dmat_t *m,
    const kalman_dmat_t *m_cov)
{
	sample_t *sample = safe_calloc(1, sizeof (*sample));

	KAL_ASSERT(vis != NULL);

	if (m != NULL)
		sample->m = kalman_dmat_copy(m);

	sample->state = kalman_get_dstate(vis->kal);
	if (KAL_IS_NULL_DMAT(*sample->state)) {
		/* Remove all samples when the Kalman filter has been reset */
		free_sample(sample);
		mutex_enter(&vis->lock);
		while ((sample = list_remove_head(&vis->samples)) != NULL)
			free_sample(sample);
		mutex_exit(&vis->lock);
		return;
	}

	mutex_enter(&vis->lock);
	while (list_count(&vis->samples) >= vis->max_samples) {
		sample_t *old_sample = list_remove_tail(&vis->samples);
		free_sample(old_sample);
	}
	list_insert_head(&vis->samples, sample);
	vis->cov = kalman_get_dcov_mat(vis->kal);
	if (m_cov != NULL) {
		free(vis->m_cov);
		vis->m_cov = kalman_dmat_copy(m_cov);
		vis->m_cov_age = 0;
	} else {
		vis->m_cov_age++;
	}
	free(vis->cont);
	vis->cont = kalman_get_dcont(vis->kal);
	mutex_exit(&vis->lock);
}

void
kalman_vis_reset(kalman_vis_t *vis)
{
	sample_t *sample;

	KAL_ASSERT(vis != NULL);

	mutex_enter(&vis->lock);
	while ((sample = list_remove_head(&vis->samples)) != NULL)
		free_sample(sample);
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
    int decimals)
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
kalman_vis_set_cov_precision(kalman_vis_t *vis, int decimals)
{
	KAL_ASSERT(vis != NULL);
	vis->cov_precision = decimals;
}

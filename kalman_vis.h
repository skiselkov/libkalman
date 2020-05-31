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

#ifndef	_LIBKALMAN_VIS_H_
#define	_LIBKALMAN_VIS_H_

#include <stdbool.h>

#include "kalman.h"

#ifdef	__cplusplus
extern "C" {
#endif

typedef struct kalman_vis_s kalman_vis_t;

kalman_vis_t *kalman_vis_alloc(kalman_t *kal, const char *name);
void kalman_vis_free(kalman_vis_t *vis);
void kalman_vis_update(kalman_vis_t *vis, const kalman_vec_t *m,
    const kalman_mat_t *m_cov);
void kalman_vis_reset(kalman_vis_t *vis);
void kalman_vis_open(kalman_vis_t *vis);
bool kalman_vis_is_open(const kalman_vis_t *vis);
void kalman_vis_set_decimals(kalman_vis_t *vis, unsigned state_var,
    unsigned decimals);
void kalman_vis_set_label(kalman_vis_t *vis, unsigned state_var,
    const char *label);
void kalman_vis_set_cov_precision(kalman_vis_t *vis, unsigned decimals);

#ifdef	__cplusplus
}
#endif

#endif	/* _LIBKALMAN_VIS_H_ */

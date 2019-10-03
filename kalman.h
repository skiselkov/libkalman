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
 * Copyright 2019 Saso Kiselkov. All rights reserved.
 */

#ifndef	_ACF_UTILS_KALMAN_H_
#define	_ACF_UTILS_KALMAN_H_

#include <acfutils/core.h>
#include <acfutils/math.h>

#ifdef	__cplusplus
extern "C" {
#endif

/*
 * WARNING: all matrices are COLUMN ORDER MAJOR.
 */

#define	KALMAN_MAX_STATE_LEN	4
#define	KALMAN_MAX_MAT_LEN	POW2(KALMAN_MAX_STATE_LEN)

typedef struct kalman_s kalman_t;

kalman_t *kalman_alloc(unsigned state_len);
void kalman_free(kalman_t *kal);

/*
 * Vectors
 */
void kalman_set_state(kalman_t *kal, const double *state);
const double *kalman_get_state(kalman_t *kal);
void kalman_set_control(kalman_t *kal, const double *control);
const double *kalman_get_control(kalman_t *kal);
void kalman_set_proc_err(kalman_t *kal, const double *proc_err);
const double *kalman_get_proc_err(kalman_t *kal);

/*
 * Matrices
 */
void kalman_set_cov_mat(kalman_t *kal, const double *cov_mat);
const double *kalman_get_cov_mat(kalman_t *kal);
void kalman_set_cov_mat_err(kalman_t *kal, const double *cont_mat_err);
const double *kalman_get_cont_mat_err(kalman_t *kal);

void kalman_set_pred_mat(kalman_t *kal, const double *pred_mat);
const double *kalman_get_pred_mat(kalman_t *kal);
void kalman_set_cont_mat(kalman_t *kal, const double *cont_mat);
const double *kalman_get_cont_mat(kalman_t *kal);

/*
 * Running the Kalman filter
 */
void kalman_step(kalman_t *kal, const double *measurement,
    const double *measurement_cov_mat);

#ifdef	__cplusplus
}
#endif

#endif	/* _ACF_UTILS_KALMAN_H_ */

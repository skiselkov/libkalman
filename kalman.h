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

#ifdef	__cplusplus
extern "C" {
#endif

/*
 * WARNING: matrices are COLUMN ORDER MAJOR. Use the KALMAN_MAT macro
 * to access individual elements of the matrix.
 * Also note, although our state is fixed at 4 elements and a 4x4 matrix,
 * this is simply for convenience of operations inside of the library. If
 * your problem requires fewer state parameters to be used, simply set
 * the corresponding vector and matrix elements to 0 to ignore them.
 */
#define	KALMAN_STATE_LEN	4
#define	KALMAN_MAT(mat, col, row)	\
	((mat)[((col) * KALMAN_STATE_LEN) + row])

typedef struct kalman_s kalman_t;
typedef double kalman_vec_t[KALMAN_STATE_LEN];
typedef double kalman_mat_t[KALMAN_STATE_LEN * KALMAN_STATE_LEN];

kalman_t *kalman_alloc(void);
void kalman_free(kalman_t *kal);

/*
 * Vectors
 */
void kalman_set_state(kalman_t *kal, const kalman_vec_t state);
void kalman_get_state(const kalman_t *kal, kalman_vec_t state);
void kalman_set_control(kalman_t *kal, const kalman_vec_t control);
void kalman_get_control(const kalman_t *kal, kalman_vec_t control);
void kalman_set_proc_err(kalman_t *kal, const kalman_vec_t proc_err);
void kalman_get_proc_err(const kalman_t *kal, kalman_vec_t proc_err);

/*
 * Matrices
 */
void kalman_set_cov_mat(kalman_t *kal, const kalman_mat_t cov_mat);
void kalman_get_cov_mat(const kalman_t *kal, kalman_mat_t cov_mat);
void kalman_set_cont_mat_err(kalman_t *kal, const kalman_mat_t cont_mat_err);
void kalman_get_cont_mat_err(const kalman_t *kal, kalman_mat_t cont_mat_err);

void kalman_set_pred_mat(kalman_t *kal, const kalman_mat_t pred_mat);
void kalman_get_pred_mat(const kalman_t *kal, kalman_mat_t pred_mat);
void kalman_set_cont_mat(kalman_t *kal, const kalman_mat_t cont_mat);
void kalman_get_cont_mat(const kalman_t *kal, kalman_mat_t cont_mat);

/*
 * Running the Kalman filter
 */
void kalman_step(kalman_t *kal, const kalman_vec_t measurement,
    const kalman_mat_t measurement_cov_mat);

#ifdef	__cplusplus
}
#endif

#endif	/* _ACF_UTILS_KALMAN_H_ */

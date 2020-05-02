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

#ifndef	_LIBKALMAN_H_
#define	_LIBKALMAN_H_

#include <stdbool.h>
#include <string.h>
#include <math.h>

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
#ifdef	KALMAN_USE_MATHC
#define	KALMAN_VEC_LEN		4
#else
#define	KALMAN_VEC_LEN		6
#endif
#define	KALMAN_VEC(...)		((kalman_vec_t){{__VA_ARGS__}})
#define	KALMAN_VECi(vec, col)	((vec).v[(col)])

typedef struct kalman_s kalman_t;

typedef struct {
	double	v[KALMAN_VEC_LEN];
} kalman_vec_t;
#define	KALMAN_ZERO_VEC		((kalman_vec_t){{ 0 }})
#define	KALMAN_NULL_VEC		((kalman_vec_t){{NAN, NAN, NAN, NAN, NAN, NAN}})
#define	KALMAN_IS_NULL_VEC(vec)	(isnan((vec).v[0]))

typedef struct {
	double	m[KALMAN_VEC_LEN * KALMAN_VEC_LEN];
} kalman_mat_t;
#define	KALMAN_ZERO_MAT	((kalman_mat_t){{ 0 }})
#ifdef	KALMAN_USE_MATHC
#define	KALMAN_NULL_MAT	((kalman_mat_t){{ \
	NAN, NAN, NAN, NAN, \
	NAN, NAN, NAN, NAN, \
	NAN, NAN, NAN, NAN, \
	NAN, NAN, NAN, NAN}})
#define	KALMAN_IDENT_MAT ((kalman_mat_t){{ \
	1, 0, 0, 0, \
	0, 1, 0, 0, \
	0, 0, 1, 0, \
	0, 0, 0, 1}})
#else	/* !defined(KALMAN_USE_MATHC) */
#define	KALMAN_NULL_MAT	((kalman_mat_t){{ \
	NAN, NAN, NAN, NAN, NAN, NAN, \
	NAN, NAN, NAN, NAN, NAN, NAN, \
	NAN, NAN, NAN, NAN, NAN, NAN, \
	NAN, NAN, NAN, NAN, NAN, NAN, \
	NAN, NAN, NAN, NAN, NAN, NAN, \
	NAN, NAN, NAN, NAN, NAN, NAN}})
#define	KALMAN_IDENT_MAT ((kalman_mat_t){{ \
	1, 0, 0, 0, 0, 0, \
	0, 1, 0, 0, 0, 0, \
	0, 0, 1, 0, 0, 0, \
	0, 0, 0, 1, 0, 0, \
	0, 0, 0, 0, 1, 0, \
	0, 0, 0, 0, 0, 1 }})
#endif	/* !defined(KALMAN_USE_MATHC) */
#define	KALMAN_IS_NULL_MAT(mat)	(isnan((mat).m[0]))

#define	KALMAN_MATxy(mat, col, row)	\
	((mat).m[((col) * KALMAN_VEC_LEN) + row])

#ifdef	KALMAN_USE_MATHC
#define	KALMAN_MAT2_BYROW(row1col1, row1col2, row2col1, row2col2) \
	((kalman_mat_t){{ \
	row1col1,	row2col1,	0,	0, \
	row1col2,	row2col2 }})
#define	KALMAN_MAT3_BYROW(\
    row1col1, row1col2, row1col3, \
    row2col1, row2col2, row2col3, \
    row3col1, row3col2, row3col3) ((kalman_mat_t){{ \
	row1col1,	row2col1,	row3col1,	0, \
	row1col2,	row2col2,	row3col2,	0, \
	row1col3,	row2col3,	row3col3}})
#define	KALMAN_MAT4_BYROW(\
    row1col1, row1col2, row1col3, row1col4, \
    row2col1, row2col2, row2col3, row2col4, \
    row3col1, row3col2, row3col3, row3col4, \
    row4col1, row4col2, row4col3, row4col4) ((kalman_mat_t){{ \
	row1col1, row2col1, row3col1, row4col1, \
	row1col2, row2col2, row3col2, row4col2, \
	row1col3, row2col3, row3col3, row4col3, \
	row1col4, row2col4, row3col4, row4col4 }})
#else	/* !defined(KALMAN_USE_MATHC) */
#define	KALMAN_MAT2_BYROW(row1col1, row1col2, row2col1, row2col2) \
	((kalman_mat_t){{ \
	row1col1,	row2col1,	0,	0,	0,	0, \
	row1col2,	row2col2 }})
#define	KALMAN_MAT3_BYROW(\
    row1col1, row1col2, row1col3, \
    row2col1, row2col2, row2col3, \
    row3col1, row3col2, row3col3) ((kalman_mat_t){{ \
	row1col1,	row2col1,	row3col1,	0,	0,	0, \
	row1col2,	row2col2,	row3col2,	0,	0,	0, \
	row1col3,	row2col3,	row3col3}})
#define	KALMAN_MAT4_BYROW(\
    row1col1, row1col2, row1col3, row1col4, \
    row2col1, row2col2, row2col3, row2col4, \
    row3col1, row3col2, row3col3, row3col4, \
    row4col1, row4col2, row4col3, row4col4) ((kalman_mat_t){{ \
	row1col1, row2col1, row3col1, row4col1, 0, 0, \
	row1col2, row2col2, row3col2, row4col2, 0, 0, \
	row1col3, row2col3, row3col3, row4col3, 0, 0, \
	row1col4, row2col4, row3col4, row4col4 }})
#define	KALMAN_MAT5_BYROW(\
    row1col1, row1col2, row1col3, row1col4, row1col5, \
    row2col1, row2col2, row2col3, row2col4, row2col5, \
    row3col1, row3col2, row3col3, row3col4, row3col5, \
    row4col1, row4col2, row4col3, row4col4, row4col5, \
    row5col1, row5col2, row5col3, row5col4, row5col5) \
    ((kalman_mat_t){{ \
	row1col1, row2col1, row3col1, row4col1, row5col1, 0, \
	row1col2, row2col2, row3col2, row4col2, row5col2, 0, \
	row1col3, row2col3, row3col3, row4col3, row5col3, 0, \
	row1col4, row2col4, row3col4, row4col4, row5col4, 0, \
	row1col5, row2col5, row3col5, row4col5, row5col5, 0}})
#define	KALMAN_MAT6_BYROW(\
    row1col1, row1col2, row1col3, row1col4, row1col5, row1col6, \
    row2col1, row2col2, row2col3, row2col4, row2col5, row2col6, \
    row3col1, row3col2, row3col3, row3col4, row3col5, row3col6, \
    row4col1, row4col2, row4col3, row4col4, row4col5, row4col6, \
    row5col1, row5col2, row5col3, row5col4, row5col5, row5col6, \
    row6col1, row6col2, row6col3, row6col4, row6col5, row6col6) \
    ((kalman_mat_t){{ \
	row1col1, row2col1, row3col1, row4col1, row5col1, row6col1, \
	row1col2, row2col2, row3col2, row4col2, row5col2, row6col2, \
	row1col3, row2col3, row3col3, row4col3, row5col3, row6col3, \
	row1col4, row2col4, row3col4, row4col4, row5col4, row6col4, \
	row1col5, row2col5, row3col5, row4col5, row5col5, row6col5, \
	row1col6, row2col6, row3col6, row4col6, row5col6, row6col6 }})
#endif	/* !defined(KALMAN_USE_MATHC) */

kalman_t *kalman_alloc(unsigned state_len);
void kalman_free(kalman_t *kal);
unsigned kalman_get_state_len(const kalman_t *kal);

/*
 * Vectors
 */
void kalman_set_state(kalman_t *kal, const kalman_vec_t *state);
kalman_vec_t kalman_get_state(const kalman_t *kal);
void kalman_set_cont(kalman_t *kal, const kalman_vec_t *control);
kalman_vec_t kalman_get_cont(const kalman_t *kal);
void kalman_set_proc_err(kalman_t *kal, const kalman_vec_t *proc_err);
kalman_vec_t kalman_get_proc_err(const kalman_t *kal);

/*
 * Matrices
 */
void kalman_set_cov_mat(kalman_t *kal, const kalman_mat_t *cov_mat);
kalman_mat_t kalman_get_cov_mat(const kalman_t *kal);
void kalman_set_cov_mat_err(kalman_t *kal, const kalman_mat_t *cov_mat_err);
kalman_mat_t kalman_get_cov_mat_err(const kalman_t *kal);

void kalman_set_pred_mat(kalman_t *kal, const kalman_mat_t *pred_mat);
kalman_mat_t kalman_get_pred_mat(const kalman_t *kal);
void kalman_set_cont_mat(kalman_t *kal, const kalman_mat_t *cont_mat);
kalman_mat_t kalman_get_cont_mat(const kalman_t *kal);

/*
 * Running the Kalman filter
 */
void kalman_step(kalman_t *kal, const kalman_vec_t *measurement,
    const kalman_mat_t *measurement_cov_mat,
    const kalman_mat_t *observation_model_p);
void kalman_step_null(kalman_t *kal);

/*
 * Utility functions
 */
void kalman_combine_measurements(double m0, double var0,
    double m1, double var1, double *m_out, double *var_out);

/*
 * Debugging the Kalman filter
 */
#define	KALMAN_PRINT_VEC(vec, state_len) \
	kalman_print_vec(#vec, (vec), (state_len))
void kalman_print_vec(const char *name, const kalman_vec_t *vec,
    unsigned state_len);
#define	KALMAN_PRINT_MAT(mat, state_len) \
	kalman_print_mat(#mat, (mat), (state_len))
void kalman_print_mat(const char *name, const kalman_mat_t *mat,
    unsigned state_len);

#define	KALMAN_DEBUG_VEC(kal, element_name, getfunc) \
	do { \
		kalman_vec_t vec = getfunc((kal)); \
		kalman_print_vec(#kal "(" element_name ")", &vec, \
		    kalman_get_state_len((kal))); \
	} while (0)
#define	KALMAN_DEBUG_MAT(kal, element_name, getfunc) \
	do { \
		kalman_mat_t mat = getfunc((kal)); \
		kalman_print_mat(#kal "(" element_name ")", &mat, \
		    kalman_get_state_len((kal))); \
	} while (0)

#define	KALMAN_DEBUG_STATE(kal) \
	KALMAN_DEBUG_VEC(kal, "state", kalman_get_state)

#define	KALMAN_DEBUG_COV_MAT(kal) \
	KALMAN_DEBUG_MAT(kal, "cov_mat", kalman_get_cov_mat)

#ifdef	__cplusplus
}
#endif

#endif	/* _LIBKALMAN_H_ */

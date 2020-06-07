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
 * Select your preferred real number type using one of the following
 * defines. If nothing is defined, libkalman defaults to double precision
 * floating point numbers.
 */
#if	defined(KALMAN_REAL_FLOAT)
typedef float kalman_real_t;
#elif	defined(KALMAN_REAL_LONG_DOUBLE)
typedef long double kalman_real_t;
#else
typedef double kalman_real_t;
#endif

/*
 * WARNING: matrices are ROW ORDER MAJOR.
 *
 * Use the KALMAN_MATyx or KALMAN_MATxy macros to access individual elements
 * of the matrix, or use the .mYX fields.
 *
 * Also note, although our state is fixed at 9 elements and a 9x9 matrix,
 * this is simply for convenience of the C-type interface of the library.
 * If your problem requires fewer state parameters to be used, simply set
 * the corresponding vector and matrix elements to 0 to ignore them. Or if
 * you need more state, set KALMAN_VEC_LEN macro from your compiler
 * options and (optionally) extend the KALMAN_IDENT_MAT macro below. This
 * macro is not used by the library, it's simply for your convenience, so
 * if you don't use it, you don't need to modify it and a simple tweak to
 * KALMAN_VEC_LEN will suffice.
 */
#ifndef	KALMAN_VEC_LEN
#define	KALMAN_VEC_LEN		9
#endif

#define	KALMAN_VEC(...)		((kalman_vec_t){{{__VA_ARGS__}}})
#define	KALMAN_VECi(vec, col)	((vec).v[(col)])

typedef struct kalman_s kalman_t;

typedef struct {
	union {
		kalman_real_t		v[KALMAN_VEC_LEN];
		struct {
			kalman_real_t	x;
			kalman_real_t	y;
			kalman_real_t	z;
			kalman_real_t	p;
			kalman_real_t	q;
			kalman_real_t	r;
			kalman_real_t	a;
			kalman_real_t	b;
			kalman_real_t	c;
		};
	};
} kalman_vec_t;
#define	KALMAN_ZERO_VEC		((kalman_vec_t){{{ 0 }}})
#define	KALMAN_NULL_VEC \
	((kalman_vec_t){{{ [0 ... (KALMAN_VEC_LEN - 1)] = NAN}}})
#define	KALMAN_IS_NULL_VEC(vec)	(isnan((vec).v[0]))

typedef struct {
	union {
		kalman_real_t		m[KALMAN_VEC_LEN * KALMAN_VEC_LEN];
		struct {
			/* first row */
			kalman_real_t	m11;
			kalman_real_t	m12;
			kalman_real_t	m13;
			kalman_real_t	m14;
			kalman_real_t	m15;
			kalman_real_t	m16;
			kalman_real_t	m17;
			kalman_real_t	m18;
			kalman_real_t	m19;
			/* second row */
			kalman_real_t	m21;
			kalman_real_t	m22;
			kalman_real_t	m23;
			kalman_real_t	m24;
			kalman_real_t	m25;
			kalman_real_t	m26;
			kalman_real_t	m27;
			kalman_real_t	m28;
			kalman_real_t	m29;
			/* third row */
			kalman_real_t	m31;
			kalman_real_t	m32;
			kalman_real_t	m33;
			kalman_real_t	m34;
			kalman_real_t	m35;
			kalman_real_t	m36;
			kalman_real_t	m37;
			kalman_real_t	m38;
			kalman_real_t	m39;
			/* fourth row */
			kalman_real_t	m41;
			kalman_real_t	m42;
			kalman_real_t	m43;
			kalman_real_t	m44;
			kalman_real_t	m45;
			kalman_real_t	m46;
			kalman_real_t	m47;
			kalman_real_t	m48;
			kalman_real_t	m49;
			/* fifth row */
			kalman_real_t	m51;
			kalman_real_t	m52;
			kalman_real_t	m53;
			kalman_real_t	m54;
			kalman_real_t	m55;
			kalman_real_t	m56;
			kalman_real_t	m57;
			kalman_real_t	m58;
			kalman_real_t	m59;
			/* sixth row */
			kalman_real_t	m61;
			kalman_real_t	m62;
			kalman_real_t	m63;
			kalman_real_t	m64;
			kalman_real_t	m65;
			kalman_real_t	m66;
			kalman_real_t	m67;
			kalman_real_t	m68;
			kalman_real_t	m69;
			/* seventh row */
			kalman_real_t	m71;
			kalman_real_t	m72;
			kalman_real_t	m73;
			kalman_real_t	m74;
			kalman_real_t	m75;
			kalman_real_t	m76;
			kalman_real_t	m77;
			kalman_real_t	m78;
			kalman_real_t	m79;
			/* eighth row */
			kalman_real_t	m81;
			kalman_real_t	m82;
			kalman_real_t	m83;
			kalman_real_t	m84;
			kalman_real_t	m85;
			kalman_real_t	m86;
			kalman_real_t	m87;
			kalman_real_t	m88;
			kalman_real_t	m89;
			/* ninth row */
			kalman_real_t	m91;
			kalman_real_t	m92;
			kalman_real_t	m93;
			kalman_real_t	m94;
			kalman_real_t	m95;
			kalman_real_t	m96;
			kalman_real_t	m97;
			kalman_real_t	m98;
			kalman_real_t	m99;
		};
	};
} kalman_mat_t;
#define	KALMAN_ZERO_MAT	((kalman_mat_t){{{ 0 }}})
#define	KALMAN_NULL_MAT	((kalman_mat_t)\
	{{{ [0 ... ((KALMAN_VEC_LEN * KALMAN_VEC_LEN) - 1)] = NAN }}})
#define	KALMAN_IDENT_MAT	((kalman_mat_t){{{ \
	1, 0, 0, 0, 0, 0, 0, 0, 0, \
	0, 1, 0, 0, 0, 0, 0, 0, 0, \
	0, 0, 1, 0, 0, 0, 0, 0, 0, \
	0, 0, 0, 1, 0, 0, 0, 0, 0, \
	0, 0, 0, 0, 1, 0, 0, 0, 0, \
	0, 0, 0, 0, 0, 1, 0, 0, 0, \
	0, 0, 0, 0, 0, 0, 1, 0, 0, \
	0, 0, 0, 0, 0, 0, 0, 1, 0, \
	0, 0, 0, 0, 0, 0, 0, 0, 1 }}})
#define	KALMAN_IS_NULL_MAT(mat)	(isnan((mat).m[0]))

#define	KALMAN_MATyx(mat, row, col)	\
	((mat).m[((row) * KALMAN_VEC_LEN) + (col)])

#define	KALMAN_MATxy(mat, col, row)	\
	((mat).m[((row) * KALMAN_VEC_LEN) + (col)])

/*
 * Allocation & destruction of the filter.
 */
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
void kalman_set_proc_noise(kalman_t *kal, const kalman_vec_t *proc_err);
kalman_vec_t kalman_get_proc_noise(const kalman_t *kal);

/*
 * Matrices
 */
void kalman_set_cov_mat(kalman_t *kal, const kalman_mat_t *cov_mat);
kalman_mat_t kalman_get_cov_mat(const kalman_t *kal);
void kalman_set_proc_noise_cov(kalman_t *kal, const kalman_mat_t *cov_mat_err);
kalman_mat_t kalman_get_proc_noise_cov(const kalman_t *kal);

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

/*
 * Utility functions
 */
void kalman_combine_s(kalman_real_t m0, kalman_real_t var0,
    kalman_real_t m1, kalman_real_t var1,
    kalman_real_t *m_out, kalman_real_t *var_out);
void kalman_combine_v(unsigned state_len,
    const kalman_vec_t *m0_in, const kalman_mat_t *cov0_in,
    const kalman_vec_t *m1_in, const kalman_mat_t *cov1_in,
    kalman_vec_t *m_out, kalman_mat_t *var_out);

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

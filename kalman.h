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

typedef struct kalman_s kalman_t;

/*
 * WARNING: matrices are ROW ORDER MAJOR.
 */

#include "kalman_deprecated.h"

typedef struct {
	unsigned	rows;
	unsigned	cols;
	kalman_real_t	dm[1];
} kalman_dmat_t;

#define	KAL_DMATyx(mat, row, col) \
	((mat).dm[(row) * (mat).cols + (col)])
#define	KAL_DMAT_SIZEOF(mat) (sizeof (kalman_dmat_t) + \
	sizeof (kalman_real_t) * ((mat).rows * (mat).cols - 1))
#define	KAL_IS_NULL_DMAT(mat)	(isnan((mat).dm[0]))

/**
 * Allocates a matrix that can be used to hold a state vector, or
 * any other matrix and passed to the Kalman filter. The matrix is
 * allocated as a single block of memory and must be freed by the
 * caller by using the standard C library free() function.
 *
 * @param rows Number of rows to be contained in the matrix
 *	(must be at least 1).
 * @param cols Number of columns to be contained in the matrix
 *	(must be at least 1).
 */
static inline kalman_dmat_t *
kalman_dmat_alloc(unsigned rows, unsigned cols)
{
	kalman_dmat_t *m = (kalman_dmat_t *)calloc(1,
	    sizeof (*m) + sizeof (kalman_real_t) * (rows * cols - 1));
	m->rows = rows;
	m->cols = cols;
	return (m);
}

/**
 * Allocates an identity matrix of `size' rows and columns. The matrix
 * is allocated as a single block of memory and must be freed by the
 * caller by using the standard C library free() function.
 */
static inline kalman_dmat_t *
kalman_dmat_alloc_ident(unsigned size)
{
	kalman_dmat_t *m = kalman_dmat_alloc(size, size);

	for (unsigned i = 0; i < size; i++)
		KAL_DMATyx(*m, i, i) = 1;

	return (m);
}

/**
 * Allocates an NAN-filled matrix of dimensions `rows' and `cols'. The
 * matrix is allocated as a single block of memory and must be freed
 * by the caller by using the standard C library free() function.
 */
static inline kalman_dmat_t *
kalman_dmat_alloc_null(unsigned rows, unsigned cols)
{
	kalman_dmat_t *m = kalman_dmat_alloc(rows, cols);

	for (unsigned row = 0; row < rows; row++) {
		for (unsigned col = 0; col < cols; col++)
			KAL_DMATyx(*m, row, col) = NAN;
	}

	return (m);
}

/**
 * Allocates a new matrix that is a copy of `mat_in'. The caller is
 * responsible for freeing the returned matrix using the standard C
 * library free() function.
 */
static inline kalman_dmat_t *
kalman_dmat_copy(const kalman_dmat_t *mat_in)
{
	kalman_dmat_t *mat_out =
	    (kalman_dmat_t *)malloc(KAL_DMAT_SIZEOF(*mat_in));
	memcpy(mat_out, mat_in, KAL_DMAT_SIZEOF(*mat_in));
	return (mat_out);
}

/*
 * Allocation & destruction of the filter.
 */
kalman_t *kalman_alloc(unsigned state_len);
void kalman_free(kalman_t *kal);
unsigned kalman_get_state_len(const kalman_t *kal);

/*
 * Manipulating the filter's current state vector and covariance matrix.
 */
void kalman_set_dstate(kalman_t *kal, const kalman_dmat_t *state);
kalman_dmat_t *kalman_get_dstate(const kalman_t *kal);

void kalman_set_dcov_mat(kalman_t *kal, const kalman_dmat_t *cov_mat);
kalman_dmat_t *kalman_get_dcov_mat(const kalman_t *kal);
/*
 * Manipulating the filter's prediction matrix.
 */
void kalman_set_dpred_mat(kalman_t *kal, const kalman_dmat_t *pred_mat);
kalman_dmat_t *kalman_get_dpred_mat(const kalman_t *kal);
/*
 * Manipulating the filter's control vector and control matrix.
 */
void kalman_set_dcont(kalman_t *kal, const kalman_dmat_t *control);
kalman_dmat_t *kalman_get_dcont(const kalman_t *kal);

void kalman_set_dcont_mat(kalman_t *kal, const kalman_dmat_t *cont_mat);
kalman_dmat_t *kalman_get_dcont_mat(const kalman_t *kal);
/*
 * Manipulating the filter's process noise vector and process noise covariance.
 */
void kalman_set_dproc_noise(kalman_t *kal, const kalman_dmat_t *proc_noise);
kalman_dmat_t *kalman_get_dproc_noise(const kalman_t *kal);

void kalman_set_dproc_noise_cov(kalman_t *kal,
    const kalman_dmat_t *proc_noise_cov);
kalman_dmat_t *kalman_get_dproc_noise_cov(const kalman_t *kal);
/*
 * Running the Kalman filter
 */
void kalman_dstep(kalman_t *kal, const kalman_dmat_t *measurement,
    const kalman_dmat_t *measurement_cov_mat,
    const kalman_dmat_t *observation_model);
/*
 * Utility functions
 */
void kalman_combine_s(kalman_real_t m0, kalman_real_t var0,
    kalman_real_t m1, kalman_real_t var1,
    kalman_real_t *m_out, kalman_real_t *var_out);

void kalman_combine(const kalman_dmat_t *m0_in, const kalman_dmat_t *cov0_in,
    const kalman_dmat_t *m1_in, const kalman_dmat_t *cov1_in,
    kalman_dmat_t *m_out, kalman_dmat_t *var_out);
/*
 * Debugging the Kalman filter
 */
#define	KALMAN_PRINT_DMAT(mat, state_len) \
	kalman_print_dmat(#mat, (mat), (state_len))
void kalman_print_dmat(const char *name, const kalman_dmat_t *mat);

#define	KALMAN_DEBUG_DMAT(kal, element_name, getfunc) \
	do { \
		kalman_dmat_t *mat = getfunc((kal)); \
		kalman_print_dmat(#kal "(" element_name ")", &mat); \
		free(mat); \
	} while (0)

#define	KALMAN_DEBUG_STATE(kal) \
	KALMAN_DEBUG_MAT(kal, "state", kalman_get_dstate)

#define	KALMAN_DEBUG_COV_MAT(kal) \
	KALMAN_DEBUG_MAT(kal, "cov_mat", kalman_get_dcov_mat)

#ifdef	__cplusplus
}
#endif

#endif	/* _LIBKALMAN_H_ */

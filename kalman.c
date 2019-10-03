/*
 * CDDL HEADER START
 *
 * This file and its contents are supplied under the terms of the
 * Common Development and Distribution License ("CDDL"), version 1.0.
 * You may only use this file in accordance with the terms of version
 * 1.0 of the CDDL.
 *
 * A full copy of the text of the CDDL should have accompanied this
 * source.  A copy of the CDDL is also available via the Internet at
 * http://www.illumos.org/license/CDDL.
 *
 * CDDL HEADER END
*/
/*
 * Copyright 2019 Saso Kiselkov. All rights reserved.
 */

#include <string.h>

#include "mathc.h"

#include "acfutils/kalman.h"
#include "acfutils/safe_alloc.h"

#define	KAL_VEC_LEN(kal) \
	(sizeof (double) * (kal)->state_len)
#define	KAL_MAT_LEN(kal) \
	(sizeof (double) * POW2((kal)->state_len))

struct kalman_s {
	unsigned	state_len;

	/* state vector */
	double		x_k[KALMAN_MAX_STATE_LEN];
	/* process error vector */
	double		u_k[KALMAN_MAX_STATE_LEN];

	/* Process covariance matrix */
	double		P_k[KALMAN_MAX_MAT_LEN];
	/* Process covariance matrix prediction error */
	double		Q_k[KALMAN_MAX_MAT_LEN];
	/* State prediction matrix */
	double		A_k[KALMAN_MAX_MAT_LEN];
	/* State prediction matrix - transposed */
	double		A_k_T[KALMAN_MAX_MAT_LEN];
	/* Control matrix */
	double		B_k[KALMAN_MAX_MAT_LEN];
};

CTASSERT(sizeof (mfloat_t) == sizeof (double));

kalman_t *
kalman_alloc(unsigned state_len)
{
	kalman_t *kal = safe_calloc(1, sizeof (*kal));

	ASSERT(state_len != 0);
	ASSERT3U(state_len, <=, KALMAN_MAX_STATE_LEN);
	ASSERT(pred_mat != NULL);

	kal->state_len = state_len;

	return (kal);
}

void
kalman_free(kalman_t *kal)
{
	free(kal);
}

void
kalman_set_state(kalman_t *kal, const double *state)
{
	ASSERT(kal != NULL);
	ASSERT(state != NULL);
	memcpy(kal->x_k, state, KAL_VEC_LEN(kal));
}

const double *
kalman_get_state(kalman_t *kal)
{
	ASSERT(kal != NULL);
	return (kal->x_k);
}

void
kalman_set_control(kalman_t *kal, const double *control)
{
	ASSERT(kal != NULL);
	ASSERT(control != NULL);
	memcpy(kal->u_k, control, KAL_VEC_LEN(kal));
}

const double *
kalman_get_control(kalman_t *kal)
{
	ASSERT(kal != NULL);
	return (kal->u_k);
}

void
kalman_set_proc_err(kalman_t *kal, const double *proc_err)
{
	ASSERT(kal != NULL);
	ASSERT(proc_err != NULL);
	memcpy(kal->u_k, proc_err, KAL_VEC_LEN(kal));
}

const double *
kalman_get_proc_err(kalman_t *kal)
{
	ASSERT(kal != NULL);
	return (kal->u_k);
}

void
kalman_set_cov_mat(kalman_t *kal, const double *cov_mat)
{
	ASSERT(kal != NULL);
	ASSERT(cov_mat != NULL);
	memcpy(kal->P_k, cov_mat, KAL_MAT_LEN(kal));
}

const double *
kalman_get_cov_mat(kalman_t *kal)
{
	ASSERT(kal != NULL);
	return (kal->P_k);
}

void
kalman_set_cov_mat_err(kalman_t *kal, const double *cov_mat_err)
{
	ASSERT(kal != NULL);
	ASSERT(cov_mat_err != NULL);
	memcpy(kal->Q_k, cov_mat_err, KAL_MAT_LEN(kal));
}

const double *
kalman_get_cov_mat_err(kalman_t *kal)
{
	ASSERT(kal != NULL);
	return (kal->Q_k);
}

void
kalman_set_pred_mat(kalman_t *kal, const double *pred_mat)
{
	ASSERT(kal != NULL);
	ASSERT(pred_mat != NULL);
	memcpy(kal->A_k, pred_mat, KAL_MAT_LEN(kal));
	mat2_transpose(kal->A_k_T, kal->A_k);
}

const double *
kalman_get_pred_mat(kalman_t *kal)
{
	ASSERT(kal != NULL);
	return (kal->A_k);
}

void
kalman_set_cont_mat(kalman_t *kal, const double *cont_mat)
{
	ASSERT(kal != NULL);
	ASSERT(cont_mat != NULL);
	memcpy(kal->B_k, cont_mat, KAL_MAT_LEN(kal));
}

const double *
kalman_get_cont_mat(kalman_t *kal)
{
	ASSERT(kal != NULL);
	return (kal->B_k);
}

void
kalman_step(kalman_t *kal, const double *measurement,
    const double *measurement_cov_mat)
{
	mfloat_t tmpvec[VEC4_SIZE], tmpvec2[VEC4_SIZE];
	mfloat_t tmpmat[MAT4_SIZE], tmpmat2[MAT4_SIZE];
	mfloat_t P_k_pred[MAT4_SIZE];
	mfloat_t x_k_pred[VEC4_SIZE];

	ASSERT(kal != NULL);
	ASSERT(measurement != NULL);
	ASSERT(measurement_cov_mat != NULL);

	/*
	 * Prediction phase:
	 *
	 * predicted state vector:
	 * ^         ^          ->
	 * x' = A  * x    + B * u
	 *  k    k    k-1    k   k
	 *
	 * predicted uncertainty:
	 *                   T
	 * P' = A  * P    * A  + Q
	 *  k    k    k-1    k    k
	 */
	vec4_multiply_mat4(tmpvec, kal->x_k, kal->A_k);
	vec4_multiply_mat4(tmpvec2, u_k, B_k);
	vec4_add(x_k_pred, tmpvec, tmpvec2);

	mat4_multiply(P_k_pred, kal->A_k, kal->P_k);
	mat4_multiply(P_k_pred, P_k_pred, kal->A_k_T);
	mat4_add(P_k_pred, P_k_pred, kal->Q_k);
	P_k_pred[1] = 0;
	P_k_pred[2] = 0;
}

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

#include <stdio.h>

#include "mathc.h"
#include "kalman.h"
#include "kalman_assert.h"

CTASSERT(sizeof (mfloat_t) == sizeof (double));

#define	MAT(mat, col, row)	KALMAN_MATxy(mat, col, row)

struct kalman_s {
	unsigned	state_len;

	/* state vector */
	kalman_vec_t	x_k;
	/* control vector */
	kalman_vec_t	u_k;
	/* process error vector */
	kalman_vec_t	w_k;

	/* Process covariance matrix */
	kalman_mat_t	P_k;
	/* Process covariance matrix prediction error */
	kalman_mat_t	Q_k;
	/* State prediction matrix */
	kalman_mat_t	A_k;
	/* State prediction matrix - transposed */
	kalman_mat_t	A_k_T;
	/* Control matrix */
	kalman_mat_t	B_k;
	/* Sensor covariance matrix */
	kalman_mat_t	R_k;
};

kalman_t *
kalman_alloc(unsigned state_len)
{
	kalman_t *kal = safe_calloc(1, sizeof (kalman_t));

	KAL_ASSERT(state_len != 0);
	KAL_ASSERT3U(state_len, <=, KALMAN_VEC_LEN);
	kal->state_len = state_len;
	/*
	 * We initialize the minimum set of vectors and matrices to null
	 * vector/matrix values, to force the user to set at least these
	 * parameters before invoking kalman_step for the first time.
	 */
	kal->x_k = KALMAN_NULL_VEC;
	kal->A_k = KALMAN_NULL_MAT;
	kal->P_k = KALMAN_NULL_MAT;

	return (kal);
}

void
kalman_free(kalman_t *kal)
{
	free(kal);
}

unsigned
kalman_get_state_len(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (kal->state_len);
}

/*
 * Sets the current state vector of the Kalman filter. Use this to set an
 * initial state.
 *
 * Theoretical background:
 * The state vector x_k represents the state of the system we want to model
 * in the Kalman filter. Each element of the vector represents one state
 * parameter. The state's parameters are updated according to the measurements
 * taken and the relationship between their covariance matrices.
 */
void
kalman_set_state(kalman_t *kal, const kalman_vec_t *state)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(state != NULL);
	kal->x_k = *state;
}

/*
 * Returns the Kalman filter's current state vector.
 */
kalman_vec_t
kalman_get_state(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (kal->x_k);
}

/*
 * Sets the Kalman filter's control vector.
 *
 * Theoretical background:
 * The control vector u_k is used to modify the prediction for the next
 * state. For example, if the Kalman filter is used to track position
 * and velocity, and we have additional knowledge of the instantaneous
 * acceleration, we can supply the acceleration component in the control
 * vector and use it to modify the predicted position and velocity at the
 * next step.
 */
void
kalman_set_cont(kalman_t *kal, const kalman_vec_t *control)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(control != NULL);
	kal->u_k = *control;
}

/*
 * Returns the Kalman filter's current control vector.
 */
kalman_vec_t
kalman_get_cont(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (kal->u_k);
}

/*
 * Sets the Kalman filter's process error vector.
 *
 * Theoretical background:
 * The process error vector w_k is a value added on top of the predicted
 * state vector to generate the predicted next state vector. It can
 * represent something like a constant error accumulated in the process
 * of predicting the next state of the system.
 */
void
kalman_set_proc_err(kalman_t *kal, const kalman_vec_t *proc_err)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(proc_err != NULL);
	kal->w_k = *proc_err;
}

/*
 * Returns the Kalman filter's process error vector.
 */
kalman_vec_t
kalman_get_proc_err(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (kal->w_k);
}

/*
 * Sets the Kalman filter's current process covariance matrix.
 *
 * Theoretical background:
 * The process covariance matrix P_k represents the uncertainty
 * (variance-covariance) in the current state of the system. You should
 * set an initial process covariance matrix to give the Kalman filter a
 * reasonable starting point. The filter will then update the covariance
 * matrix every step, depending on the covariance matrices of the
 * prediction step and measurement.
 */
void
kalman_set_cov_mat(kalman_t *kal, const kalman_mat_t *cov_mat)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(cov_mat != NULL);
	kal->P_k = *cov_mat;
}

/*
 * Returns the Kalman filter's current process covariance matrix.
 */
kalman_mat_t
kalman_get_cov_mat(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (kal->P_k);
}

/*
 * Sets the Kalman filter's current covariance matrix prediction error.
 *
 * Theoretical background:
 * In the prediction step, the Kalman filter predicts the next process
 * covariance matrix (i.e. "process uncertainty"). During this stage,
 * we may want to add extra uncertainty to the predictor. For example,
 * if we want to keep the predictor from reaching an eternal "low uncertainty"
 * from its own predictions, thus excessively relying on predictions, instead
 * of measurements, we can add extra uncertainty each prediction step.
 * In essence, this means the predictor is becoming more uncertain about
 * the real state of the system if it keeps on relying on the predictions
 * a lot. This in turn will mean that measurements will, by comparison, start
 * having more of an effect, thus resulting in faster response times to
 * changes in external conditions.
 */
void
kalman_set_cov_mat_err(kalman_t *kal, const kalman_mat_t *cov_mat_err)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(cov_mat_err != NULL);
	kal->Q_k = *cov_mat_err;
}

/*
 * Returns the Kalman filter's process covariance covariance matrix
 * prediction error.
 */
kalman_mat_t
kalman_get_cov_mat_err(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (kal->Q_k);
}

/*
 * Sets the Kalman filter's prediction matrix.
 *
 * Theoretical background:
 * This matrix is used to evolve the state of the system to the next
 * step. It is used to multiply the current state vector to arrive at
 * the next predicted system state.
 */
void
kalman_set_pred_mat(kalman_t *kal, const kalman_mat_t *pred_mat)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(pred_mat != NULL);
	kal->A_k = *pred_mat;
	mat4_transpose(kal->A_k_T.m, kal->A_k.m);
}

/*
 * Returns the Kalman filter's current prediction matrix.
 */
kalman_mat_t
kalman_get_pred_mat(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (kal->A_k);
}

/*
 * Sets the Kalman filter's control matrix.
 *
 * Theoretical background:
 * The control matrix adapts the control vector and transforms the
 * values in the control vector to apply properly to the corresponding
 * state vector components. For example, in the case of a position-velocity
 * Kalman filter with acceleration in the control vector, the control matrix
 * adds the acceleration-induced changes to the position and velocity
 * state parameters via the standard "d=1/2at^2" and "v=at" equations.
 * So the acceleration would be contained in the control vector, and the
 * control matrix would contain "1/2t^2" and "t" in the diagonal elements
 * of the matrix.
 */
void
kalman_set_cont_mat(kalman_t *kal, const kalman_mat_t *cont_mat)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(cont_mat != NULL);
	kal->B_k = *cont_mat;
}

/*
 * Returns the Kalman filter's control matrix.
 */
kalman_mat_t
kalman_get_cont_mat(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (kal->B_k);
}

static inline void
mat4_add(mfloat_t outmat4[MAT4_SIZE], const mfloat_t m1[MAT4_SIZE],
    const mfloat_t m2[MAT4_SIZE])
{
	for (int i = 0; i < MAT4_SIZE; i++)
		outmat4[i] = m1[i] + m2[i];
}

static void
mat4_subtract(mfloat_t outmat4[MAT4_SIZE], const mfloat_t m1[MAT4_SIZE],
    const mfloat_t m2[MAT4_SIZE])
{
	for (int i = 0; i < MAT4_SIZE; i++)
		outmat4[i] = m1[i] - m2[i];
}

/*
 * Performs one step of the Kalman filter, updating the filter's internal
 * state vector and process covariance matrix with a new measurement and
 * measurement covariance matrix.
 *
 * Theoretical background:
 * The Kalman filter combines the current state with a prediction and
 * measurement to arrive at the next state of the system. The state's
 * preference for either using its internal predicted state, or the
 * measured state depends on the ratio of the predicted process
 * covariance matrix (calculated internally by the filter) and the
 * measurement covariance matrix supplied in this call.
 */
void
kalman_step(kalman_t *kal, const kalman_vec_t *measurement,
    const kalman_mat_t *measurement_cov_mat,
    const kalman_mat_t *observation_model_p)
{
	kalman_vec_t x_k_pred, tmpvec, tmpvec2;
	kalman_mat_t P_k_pred, tmpmat, tmpmat2;
	kalman_mat_t K_numer, K_denom, K_denom_inv;
	kalman_mat_t K = KALMAN_ZERO_MAT;
	kalman_mat_t observation_model, observation_model_T;

	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(!KALMAN_IS_NULL_VEC(kal->x_k));
	KAL_ASSERT(!KALMAN_IS_NULL_MAT(kal->A_k));
	KAL_ASSERT(!KALMAN_IS_NULL_MAT(kal->P_k));
	KAL_ASSERT(measurement != NULL);
	KAL_ASSERT(measurement_cov_mat != NULL);

	if (observation_model_p != NULL)
		observation_model = *observation_model_p;
	else
		observation_model = KALMAN_IDENT_MAT;
	mat4_transpose(observation_model_T.m, observation_model.m);
	/*
	 * Prediction phase:
	 *
	 * predicted state vector:
	 *                      ->
	 * x' = A  * x    + B * u  + w
	 *  k    k    k-1    k   k    k
	 *
	 * predicted uncertainty:
	 *                   T
	 * P' = A  * P    * A  + Q
	 *  k    k    k-1    k    k
	 */
	vec4_multiply_mat4(tmpvec.v, kal->x_k.v, kal->A_k.m);
	vec4_multiply_mat4(tmpvec2.v, kal->u_k.v, kal->B_k.m);
	vec4_add(x_k_pred.v, tmpvec.v, tmpvec2.v);
	vec4_add(x_k_pred.v, x_k_pred.v, kal->w_k.v);

	mat4_multiply(tmpmat.m, kal->A_k.m, kal->P_k.m);
	mat4_multiply(tmpmat2.m, tmpmat.m, kal->A_k_T.m);
	mat4_add(P_k_pred.m, tmpmat2.m, kal->Q_k.m);

	/*
	 * Compute the Kalman gain:
	 *               T
	 *         P' * H
	 *          k    k
	 * K = -----------------
	 *                T
	 *     H  * P' * H  + R
	 *      k    k    k    k
	 */
	/* numerator portion */
	mat4_multiply(K_numer.m, P_k_pred.m, observation_model_T.m);
	/* denominator portion */
	mat4_multiply(tmpmat.m, observation_model.m, K_numer.m);
	mat4_add(K_denom.m, tmpmat.m, measurement_cov_mat->m);
	/*
	 * Combine the numerator & denominator. Careful about inversions,
	 * we need to set any elements outside of the state size to unity,
	 * otherwise the inversion will result in div-by-zero.
	 */
	switch (kal->state_len) {
	case 1:
		KALMAN_MATxy(K_denom, 1, 1) = 1;
		KALMAN_MATxy(K_denom, 2, 2) = 1;
		KALMAN_MATxy(K_denom, 3, 3) = 1;
		break;
	case 2:
		KALMAN_MATxy(K_denom, 2, 2) = 1;
		KALMAN_MATxy(K_denom, 3, 3) = 1;
		break;
	case 3:
		KALMAN_MATxy(K_denom, 3, 3) = 1;
		break;
	case 4:
		break;
	default:
		KAL_VERIFY(0);
	}
	mat4_inverse(K_denom_inv.m, K_denom.m);
	mat4_multiply(K.m, K_numer.m, K_denom_inv.m);

	/*
	 * Update the estimate via the measurement:
	 *                /            \
	 * x  = x' + K * ( z  - H * x'  )
	 *  k    k        \ k    k   k /
	 */
	vec4_multiply_mat4(tmpvec.v, x_k_pred.v, observation_model.m);
	vec4_subtract(tmpvec2.v, measurement->v, tmpvec.v);
	vec4_multiply_mat4(tmpvec.v, tmpvec2.v, K.m);
	vec4_add(kal->x_k.v, x_k_pred.v, tmpvec.v);
	/*
	 * P  = P' - K * H  * P'
	 *  k    k        k    k
	 */
	mat4_multiply(tmpmat.m, K.m, observation_model.m);
	mat4_multiply(tmpmat2.m, tmpmat.m, P_k_pred.m);
	mat4_subtract(kal->P_k.m, P_k_pred.m, tmpmat2.m);

	KAL_ASSERT3F(MAT(kal->P_k, 0, 0), >=, 0);
	KAL_ASSERT(isfinite(KALMAN_VECi(kal->x_k, 0)));
}

/*
 * This function is the null-update version of kalman_step. It can be used
 * for cases where we want to step the filter's prediction forward by some
 * fixed time quantum, but we don't have a new measurement to incorporate.
 * So we simply evolve the filter forward based on the prediction parameters.
 * This is effectively the same as kalman_step with an infinite covariance
 * (i.e. measured value is completely unreliable and thus only the
 * prediction is used).
 */
void
kalman_step_null(kalman_t *kal)
{
	kalman_vec_t tmpvec, tmpvec2, x_k;
	kalman_mat_t tmpmat, tmpmat2;

	KAL_ASSERT(kal != NULL);
	/*
	 * In the null update model, our predicition simply becomes the new
	 * state of the filter. This is because a null update behaves as if
	 * the measurement covariance is infinite, resulting in the Kalman
	 * gain K=0.
	 *
	 * new state vector:
	 *                      ->
	 * x  = A  * x    + B * u  + w
	 *  k    k    k-1    k   k    k
	 *
	 * new uncertainty:
	 *                   T
	 * P  = A  * P    * A  + Q
	 *  k    k    k-1    k    k
	 */
	vec4_multiply_mat4(tmpvec.v, kal->x_k.v, kal->A_k.m);
	vec4_multiply_mat4(tmpvec2.v, kal->u_k.v, kal->B_k.m);
	vec4_add(x_k.v, tmpvec.v, tmpvec2.v);
	vec4_add(kal->x_k.v, x_k.v, kal->w_k.v);

	mat4_multiply(tmpmat.m, kal->A_k.m, kal->P_k.m);
	mat4_multiply(tmpmat2.m, tmpmat.m, kal->A_k_T.m);
	mat4_add(kal->P_k.m, tmpmat2.m, kal->Q_k.m);
}

void
kalman_print_vec(const char *name, const kalman_vec_t *vec, unsigned state_len)
{
	printf("%s = (", name);
	for (unsigned i = 0; i < state_len; i++) {
		printf("%f%s", KALMAN_VECi(*vec, i),
		    (i + 1 < state_len) ? "\t": "");
	}
	printf(")\n");
}

void
kalman_print_mat(const char *name, const kalman_mat_t *mat, unsigned state_len)
{
	int len = strlen(name) + 4;
	char leadspace[len + 1];

	printf("%s = (", name);
	memset(leadspace, ' ', len);
	leadspace[len] = '\0';
	leadspace[len - 1] = '(';
	for (unsigned r = 0; r < state_len; r++) {
		if (r > 0)
			printf("%s", leadspace);
		for (unsigned c = 0; c < state_len; c++) {
			printf("%f%s", KALMAN_MATxy(*mat, c, r),
			    (c + 1 < state_len) ? "\t": "");
		}
		printf(")\n");
	}
}

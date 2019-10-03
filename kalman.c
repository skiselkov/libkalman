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

#include <acfutils/safe_alloc.h>

#include "mathc.h"
#include "kalman.h"

CTASSERT(sizeof (mfloat_t) == sizeof (double));

#define	MAT(mat, col, row)	KALMAN_MAT(mat, col, row)

struct kalman_s {
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
};

kalman_t *
kalman_alloc(void)
{
	return (safe_calloc(1, sizeof (kalman_t)));
}

void
kalman_free(kalman_t *kal)
{
	free(kal);
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
kalman_set_state(kalman_t *kal, const kalman_vec_t state)
{
	ASSERT(kal != NULL);
	ASSERT(state != NULL);
	memcpy(kal->x_k, state, sizeof (kalman_vec_t));
}

/*
 * Returns the Kalman filter's current state vector.
 */
void
kalman_get_state(const kalman_t *kal, kalman_vec_t state)
{
	ASSERT(kal != NULL);
	ASSERT(state != NULL);
	memcpy(state, kal->x_k, sizeof (kalman_vec_t));
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
kalman_set_control(kalman_t *kal, const kalman_vec_t control)
{
	ASSERT(kal != NULL);
	ASSERT(control != NULL);
	memcpy(kal->u_k, control, sizeof (kalman_vec_t));
}

/*
 * Returns the Kalman filter's current control vector.
 */
void
kalman_get_control(const kalman_t *kal, kalman_vec_t control)
{
	ASSERT(kal != NULL);
	ASSERT(control != NULL);
	memcpy(control, kal->u_k, sizeof (kalman_vec_t));
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
kalman_set_proc_err(kalman_t *kal, const kalman_vec_t proc_err)
{
	ASSERT(kal != NULL);
	ASSERT(proc_err != NULL);
	memcpy(kal->w_k, proc_err, sizeof (kalman_vec_t));
}

/*
 * Returns the Kalman filter's process error vector.
 */
void
kalman_get_proc_err(const kalman_t *kal, kalman_vec_t proc_err)
{
	ASSERT(kal != NULL);
	ASSERT(proc_err != NULL);
	memcpy(proc_err, kal->w_k, sizeof (kalman_vec_t));
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
kalman_set_cov_mat(kalman_t *kal, const kalman_mat_t cov_mat)
{
	ASSERT(kal != NULL);
	ASSERT(cov_mat != NULL);
	memcpy(kal->P_k, cov_mat, sizeof (kalman_mat_t));
}

/*
 * Returns the Kalman filter's current process covariance matrix.
 */
void
kalman_get_cov_mat(const kalman_t *kal, kalman_mat_t cov_mat)
{
	ASSERT(kal != NULL);
	ASSERT(cov_mat != NULL);
	memcpy(cov_mat, kal->P_k, sizeof (kalman_mat_t));
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
kalman_set_cov_mat_err(kalman_t *kal, const kalman_mat_t cov_mat_err)
{
	ASSERT(kal != NULL);
	ASSERT(cov_mat_err != NULL);
	memcpy(kal->Q_k, cov_mat_err, sizeof (kalman_mat_t));
}

/*
 * Returns the Kalman filter's process covariance covariance matrix
 * prediction error.
 */
void
kalman_get_cov_mat_err(const kalman_t *kal, kalman_mat_t cov_mat_err)
{
	ASSERT(kal != NULL);
	ASSERT(cov_mat_err != NULL);
	memcpy(cov_mat_err, kal->Q_k, sizeof (kalman_mat_t));
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
kalman_set_pred_mat(kalman_t *kal, const kalman_mat_t pred_mat)
{
	ASSERT(kal != NULL);
	ASSERT(pred_mat != NULL);
	memcpy(kal->A_k, pred_mat, sizeof (kalman_mat_t));
	mat2_transpose(kal->A_k_T, kal->A_k);
}

/*
 * Returns the Kalman filter's current prediction matrix.
 */
void
kalman_get_pred_mat(const kalman_t *kal, kalman_mat_t pred_mat)
{
	ASSERT(kal != NULL);
	ASSERT(pred_mat != NULL);
	memcpy(pred_mat, kal->A_k, sizeof (kalman_mat_t));
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
kalman_set_cont_mat(kalman_t *kal, const kalman_mat_t cont_mat)
{
	ASSERT(kal != NULL);
	ASSERT(cont_mat != NULL);
	memcpy(kal->B_k, cont_mat, sizeof (kalman_mat_t));
}

/*
 * Returns the Kalman filter's control matrix.
 */
void
kalman_get_cont_mat(const kalman_t *kal, kalman_mat_t cont_mat)
{
	ASSERT(kal != NULL);
	ASSERT(cont_mat != NULL);
	memcpy(cont_mat, kal->B_k, sizeof (kalman_mat_t));
}

static inline void
mat4_add(mfloat_t outmat2[MAT4_SIZE], const mfloat_t m1[MAT4_SIZE],
    const mfloat_t m2[MAT4_SIZE])
{
	for (int i = 0; i < MAT4_SIZE; i++)
		outmat2[i] = m1[i] + m2[i];
}

static void
mat4_subtract(mfloat_t outmat2[MAT4_SIZE], const mfloat_t m1[MAT4_SIZE],
    const mfloat_t m2[MAT4_SIZE])
{
	for (int i = 0; i < MAT4_SIZE; i++)
		outmat2[i] = m1[i] - m2[i];
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
kalman_step(kalman_t *kal, const kalman_vec_t measurement,
    const kalman_mat_t measurement_cov_mat)
{
	kalman_vec_t x_k_pred, tmpvec, tmpvec2;
	kalman_mat_t P_k_pred, K, tmpmat;

	ASSERT(kal != NULL);
	ASSERT(measurement != NULL);
	ASSERT(measurement_cov_mat != NULL);

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
	vec4_multiply_mat4(tmpvec, kal->x_k, kal->A_k);
	vec4_multiply_mat4(tmpvec2, kal->u_k, kal->B_k);
	vec4_add(x_k_pred, tmpvec, tmpvec2);
	vec4_add(x_k_pred, x_k_pred, kal->w_k);

	mat4_multiply(P_k_pred, kal->A_k, kal->P_k);
	mat4_multiply(P_k_pred, P_k_pred, kal->A_k_T);
	mat4_add(P_k_pred, P_k_pred, kal->Q_k);

	/* Set all off-diagonal elements to zero */
	MAT(P_k_pred, 0, 1) = 0;
	MAT(P_k_pred, 0, 2) = 0;
	MAT(P_k_pred, 0, 3) = 0;
	MAT(P_k_pred, 1, 0) = 0;
	MAT(P_k_pred, 1, 2) = 0;
	MAT(P_k_pred, 1, 3) = 0;
	MAT(P_k_pred, 2, 0) = 0;
	MAT(P_k_pred, 2, 1) = 0;
	MAT(P_k_pred, 2, 3) = 0;
	MAT(P_k_pred, 3, 0) = 0;
	MAT(P_k_pred, 3, 1) = 0;
	MAT(P_k_pred, 3, 2) = 0;

	/*
	 * Compute the Kalman gain:
	 *               T
	 *         P' * H
	 *          k    k
	 * K = ----------------
	 *                T
	 *     H  * P' * H  + R
	 *      k    k    k
	 */
	mat4_add(tmpmat, P_k_pred, measurement_cov_mat);

	memset(K, 0, sizeof (K));
	MAT(K, 0, 0) = MAT(P_k_pred, 0, 0) / MAT(tmpmat, 0, 0);
	MAT(K, 1, 1) = MAT(P_k_pred, 1, 1) / MAT(tmpmat, 1, 1);
	MAT(K, 2, 2) = MAT(P_k_pred, 2, 2) / MAT(tmpmat, 2, 2);
	MAT(K, 3, 3) = MAT(P_k_pred, 3, 3) / MAT(tmpmat, 3, 3);

	/*
	 * Update the estimate via the measurement:
	 *
	 *                /            \
	 * x  = x' + K * ( z  - H * x'  )
	 *  k    k        \ k    k   k /
	 *
	 */
	vec4_subtract(tmpvec, (mfloat_t *)measurement, x_k_pred);
	vec4_multiply_mat4(tmpvec, tmpvec, K);
	vec4_add(kal->x_k, x_k_pred, tmpvec);
	/*
	 * P  = P' - K * H  * P'
	 *  k    k        k    k
	 */
	mat4_multiply(tmpmat, K, P_k_pred);
	mat4_subtract(kal->P_k, P_k_pred, tmpmat);

	ASSERT3F(kal->P_k[0], >=, 0);
	ASSERT(isfinite(kal->x_k[0]));
}

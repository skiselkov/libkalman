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
 * Copyright 2020 Saso Kiselkov. All rights reserved.
 */

#include <stdio.h>

#include "kalman_assert.h"
#include "kalman.h"

#include "Eigen/Dense"
using Eigen::MatrixXd;
using Eigen::VectorXd;

#define	KAL_MAT_COPYIN(eig_mat, kal_mat, dim) \
	do { \
		const kalman_mat_t *k_mat = (kal_mat); \
		for (unsigned col = 0; col < (dim); col++) { \
			for (unsigned row = 0; row < (dim); row++) { \
				(eig_mat)(row, col) = \
				    KALMAN_MATxy(*k_mat, col, row); \
			} \
		} \
	} while (0)

#define	KAL_MAT_COPYOUT(kal_mat, eig_mat, dim) \
	do { \
		kalman_mat_t *k_mat = (kal_mat); \
		memset(k_mat, 0, sizeof (*k_mat)); \
		for (unsigned col = 0; col < dim; col++) { \
			for (unsigned row = 0; row < dim; row++) { \
				KALMAN_MATxy(*k_mat, row, col) = \
				    (eig_mat)(row, col); \
			} \
		} \
	} while (0)

#define	KAL_VEC_COPYIN(eig_vec, kal_vec, dim) \
	do { \
		const kalman_vec_t *k_vec = (kal_vec); \
		for (unsigned row = 0; row < dim; row++) \
			(eig_vec)(row) = k_vec->v[row]; \
	} while (0)

#define	KAL_VEC_COPYOUT(kal_vec, eig_vec, dim) \
	do { \
		kalman_vec_t *k_vec = (kal_vec); \
		memset(k_vec, 0, sizeof (*k_vec)); \
		for (unsigned row = 0; row < dim; row++) \
			k_vec->v[row] = (eig_vec)(row); \
	} while (0)

#define	MAT_INIT_VALUE(eig_mat, dim, value) \
	do { \
		for (unsigned row = 0; row < dim; row++) { \
			for (unsigned col = 0; col < dim; col++) \
				(eig_mat)(row, col) = (value); \
		} \
	} while (0)

#define	VEC_INIT_VALUE(eig_vec, dim, value) \
	do { \
		for (unsigned row = 0; row < dim; row++) \
			(eig_vec)(row) = (value); \
	} while (0)

struct kalman_s {
	unsigned	state_len;

	VectorXd	x_k;	/* state vector */
	VectorXd	u_k;	/* control vector */
	VectorXd	w_k;	/* process error vector */

	MatrixXd	P_k;	/* Process covariance matrix */
	MatrixXd	Q_k;	/* Process cov. matrix prediction error */
	MatrixXd	A_k;	/* State prediction matrix */
	MatrixXd	B_k;	/* Control matrix */
};

kalman_t *
kalman_alloc(unsigned state_len)
{
	kalman_t *kal = new kalman_t;

	KAL_ASSERT(state_len != 0);
	KAL_ASSERT3U(state_len, <=, KALMAN_VEC_LEN);
	memset(kal, 0, sizeof (*kal));
	kal->state_len = state_len;
	/*
	 * We init the state vector, covariance matrix and prediction matrix
	 * to NAN values. This forces the caller to set them before calling
	 * kalman_step for the first time.
	 */
	kal->x_k = VectorXd(state_len);
	VEC_INIT_VALUE(kal->x_k, state_len, NAN);
	kal->u_k = VectorXd::Zero(state_len);
	kal->w_k = VectorXd::Zero(state_len);

	kal->P_k = MatrixXd(state_len, state_len);
	MAT_INIT_VALUE(kal->P_k, state_len, NAN);
	kal->Q_k = MatrixXd::Zero(state_len, state_len);
	kal->A_k = MatrixXd(state_len, state_len);
	MAT_INIT_VALUE(kal->A_k, state_len, NAN);
	kal->B_k = MatrixXd::Zero(state_len, state_len);

	return (kal);
}

void
kalman_free(kalman_t *kal)
{
	delete kal;
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
	KAL_VEC_COPYIN(kal->x_k, state, kal->state_len);
}

/*
 * Returns the Kalman filter's current state vector.
 */
kalman_vec_t
kalman_get_state(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	kalman_vec_t v;
	KAL_VEC_COPYOUT(&v, kal->x_k, kal->state_len);
	return (v);
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
	KAL_VEC_COPYIN(kal->u_k, control, kal->state_len);
}

/*
 * Returns the Kalman filter's current control vector.
 */
kalman_vec_t
kalman_get_cont(const kalman_t *kal)
{
	kalman_vec_t v;
	KAL_ASSERT(kal != NULL);
	KAL_VEC_COPYOUT(&v, kal->u_k, kal->state_len);
	return (v);
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
	KAL_VEC_COPYIN(kal->w_k, proc_err, kal->state_len);
}

/*
 * Returns the Kalman filter's process error vector.
 */
kalman_vec_t
kalman_get_proc_err(const kalman_t *kal)
{
	kalman_vec_t v;
	KAL_ASSERT(kal != NULL);
	KAL_VEC_COPYOUT(&v, kal->w_k, kal->state_len);
	return (v);
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
	KAL_MAT_COPYIN(kal->P_k, cov_mat, kal->state_len);
}

/*
 * Returns the Kalman filter's current process covariance matrix.
 */
kalman_mat_t
kalman_get_cov_mat(const kalman_t *kal)
{
	kalman_mat_t m;
	KAL_ASSERT(kal != NULL);
	KAL_MAT_COPYOUT(&m, kal->P_k, kal->state_len);
	return (m);
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
	KAL_MAT_COPYIN(kal->Q_k, cov_mat_err, kal->state_len);
}

/*
 * Returns the Kalman filter's process covariance covariance matrix
 * prediction error.
 */
kalman_mat_t
kalman_get_cov_mat_err(const kalman_t *kal)
{
	kalman_mat_t m;
	KAL_ASSERT(kal != NULL);
	KAL_MAT_COPYOUT(&m, kal->Q_k, kal->state_len);
	return (m);
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
	KAL_MAT_COPYIN(kal->A_k, pred_mat, kal->state_len);
}

/*
 * Returns the Kalman filter's current prediction matrix.
 */
kalman_mat_t
kalman_get_pred_mat(const kalman_t *kal)
{
	kalman_mat_t m;
	KAL_ASSERT(kal != NULL);
	KAL_MAT_COPYOUT(&m, kal->A_k, kal->state_len);
	return (m);
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
	KAL_MAT_COPYIN(kal->B_k, cont_mat, kal->state_len);
}

/*
 * Returns the Kalman filter's control matrix.
 */
kalman_mat_t
kalman_get_cont_mat(const kalman_t *kal)
{
	kalman_mat_t m;
	KAL_ASSERT(kal != NULL);
	KAL_MAT_COPYOUT(&m, kal->B_k, kal->state_len);
	return (m);
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
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(!isnan(kal->x_k(0)));
	KAL_ASSERT(!isnan(kal->P_k(0, 0)));
	KAL_ASSERT(!isnan(kal->A_k(0, 0)));

	VectorXd m(kal->state_len);
	MatrixXd m_cov_mat(kal->state_len, kal->state_len);
	MatrixXd obsv_model(kal->state_len, kal->state_len);
	MatrixXd obsv_model_T(kal->state_len, kal->state_len);
	VectorXd x_k_pred(kal->state_len);
	MatrixXd P_k_pred(kal->state_len, kal->state_len);
	MatrixXd tmpmat(kal->state_len, kal->state_len);
	MatrixXd K(kal->state_len, kal->state_len);

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
	x_k_pred = (kal->A_k * kal->x_k) + (kal->B_k * kal->u_k) + kal->w_k;
	P_k_pred = (kal->A_k * kal->P_k * kal->A_k.transpose()) + kal->Q_k;

	if (measurement == NULL) {
		KAL_ASSERT(measurement_cov_mat == NULL);
		kal->x_k = x_k_pred;
		kal->P_k = P_k_pred;
		return;
	}
	KAL_ASSERT(!KALMAN_IS_NULL_VEC(*measurement));
	KAL_ASSERT(!KALMAN_IS_NULL_MAT(*measurement_cov_mat));
	KAL_ASSERT(measurement_cov_mat != NULL);

	KAL_VEC_COPYIN(m, measurement, kal->state_len);
	KAL_MAT_COPYIN(m_cov_mat, measurement_cov_mat, kal->state_len);

	if (observation_model_p != NULL)
		KAL_MAT_COPYIN(obsv_model, observation_model_p, kal->state_len);
	else
		obsv_model.setIdentity();
	obsv_model_T = obsv_model.transpose();
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
	tmpmat = P_k_pred * obsv_model_T;
	K = tmpmat * (obsv_model * tmpmat + m_cov_mat).inverse();
	/*
	 * Update the estimate via the measurement:
	 *                /            \
	 * x  = x' + K * ( z  - H * x'  )
	 *  k    k        \ k    k   k /
	 */
	kal->x_k = x_k_pred + K * (m - obsv_model * x_k_pred);
	/*
	 * P  = P' - K * H  * P'
	 *  k    k        k    k
	 */
	kal->P_k = P_k_pred - K * obsv_model * P_k_pred;

	KAL_ASSERT3F(kal->P_k(0, 0), >=, 0);
	KAL_ASSERT(isfinite(kal->x_k(0)));
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

/*
 * Takes two scalar measurements and combines them together according to
 * their variances (sigma^2), producing a new combined measurement and
 * variance. Use this if you want to update your Kalman filter with
 * multiple measurements in a single step.
 *
 * @param m0 First measurement.
 * @param var0 Variance of first measurement.
 * @param m1 Second measurement.
 * @param var1 Variance of second measurement.
 * @param m_out Output that will be populated with the combined measurement.
 * @param var_out Output that will be populated with the combined variance.
 */
void
kalman_combine_s(double m0, double var0, double m1, double var1,
    double *m_out, double *var_out)
{
	double k;

	KAL_ASSERT(m_out != NULL);
	KAL_ASSERT(var_out != NULL);
	/*
	 * Two normally distributed measurements m0 and m1 with variances
	 * var0 and var1 can be combined using the following equations to
	 * construct a combined measurement m' and var':
	 *
	 *           var0 * (m1 - m0)
	 * m' = m0 + ----------------
	 *             var0 + var1
	 *
	 *                 var0^2
	 * var' = var0 - -----------
	 *               var0 + var1
	 *
	 * Simplifying:
	 *
	 * k = var0 / (var0 + var1)
	 * m' = m0 + k(m1 - m0)
	 * var' = var0 - k * var0
	 */
	k = var0 / (var0 + var1);
	*m_out = m0 + k * (m1 - m0);
	*var_out = var0 - k * var1;
}

/*
 * Same as kalman_combine_s, but uses vectors & matrices for the measurements
 * and covariances.
 *
 * @param state_len State length of the measurements. This is the length of
 *	the vectors and dimension of the covariance matrices.
 * @param m0_in First measurement.
 * @param cov0_in Covariance of first measurement.
 * @param m1_in Second measurement.
 * @param cov1_in Covariance of second measurement.
 * @param m_out Output that will be populated with the combined measurement.
 * @param cov_out Output that will be populated with the combined covariance.
 */
void kalman_combine_v(unsigned state_len,
    const kalman_vec_t *m0_in, const kalman_mat_t *cov0_in,
    const kalman_vec_t *m1_in, const kalman_mat_t *cov1_in,
    kalman_vec_t *m_out, kalman_mat_t *cov_out)
{
	KAL_ASSERT(state_len != 0);
	KAL_ASSERT(state_len <= KALMAN_VEC_LEN);

	VectorXd m0(state_len), m1(state_len), m(state_len);
	MatrixXd cov0(state_len, state_len), cov1(state_len, state_len);
	MatrixXd cov(state_len, state_len);
	MatrixXd K(state_len, state_len);

	KAL_ASSERT(m0_in != NULL);
	KAL_ASSERT(!KALMAN_IS_NULL_VEC(*m0_in));
	KAL_ASSERT(cov0_in != NULL);
	KAL_ASSERT(!KALMAN_IS_NULL_MAT(*cov0_in));
	KAL_ASSERT(m1_in != NULL);
	KAL_ASSERT(!KALMAN_IS_NULL_VEC(*m1_in));
	KAL_ASSERT(cov1_in != NULL);
	KAL_ASSERT(!KALMAN_IS_NULL_MAT(*cov1_in));
	KAL_ASSERT(m_out != NULL);
	KAL_ASSERT(cov_out != NULL);

	KAL_VEC_COPYIN(m0, m0_in, state_len);
	KAL_MAT_COPYIN(cov0, cov0_in, state_len);
	KAL_VEC_COPYIN(m1, m1_in, state_len);
	KAL_MAT_COPYIN(cov1, cov1_in, state_len);

	K = cov0 * (cov0 + cov1).inverse();
	/*
	 * m' = m0 + K(m1 - m0)
	 */
	m = m0 + K * (m1 - m0);
	KAL_VEC_COPYOUT(m_out, m, state_len);
	KAL_ASSERT(!KALMAN_IS_NULL_VEC(*m_out));
	/*
	 * var' = var0 - K * var0
	 */
	cov = cov0 - K * cov0;
	KAL_MAT_COPYOUT(cov_out, cov, state_len);
	KAL_ASSERT(!KALMAN_IS_NULL_MAT(*cov_out));
}

void
kalman_print_vec(const char *name, const kalman_vec_t *vec, unsigned state_len)
{
	printf("%s = (", name);
	for (unsigned i = 0; i < state_len; i++) {
		printf("%-11.5f%s", KALMAN_VECi(*vec, i),
		    (i + 1 < state_len) ? " " : "");
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
			printf("%-11.5f%s", KALMAN_MATxy(*mat, c, r),
			    (c + 1 < state_len) ? " " : "");
		}
		printf(")\n");
	}
}

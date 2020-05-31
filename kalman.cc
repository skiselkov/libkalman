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

#include "Eigen/Dense"

#include "kalman_assert.h"
#include "kalman.h"

typedef Eigen::Matrix<kalman_real_t, Eigen::Dynamic, Eigen::Dynamic> KalMat;
typedef Eigen::Matrix<kalman_real_t, Eigen::Dynamic, 1> KalVec;

#ifdef	KAL_HEAVY_DEBUG
#define	CHECK_MAT(mat) \
	do { \
		for (int col = 0, cols = (mat).cols(); col < cols; col++) { \
			for (int row = 0, rows = (mat).rows(); row < rows; \
			    row++) {\
				KAL_ASSERT(!isnan((mat)(row, col))); \
			} \
		} \
	} while (0)
#else	/* !defined(KAL_HEAVY_DEBUG) */
#define	CHECK_MAT(mat)
#endif	/* !defined(KAL_HEAVY_DEBUG) */

struct kalman_s {
	unsigned	state_len;

	KalVec		x_k;	/* state vector */
	KalVec		u_k;	/* control vector */
	KalVec		w_k;	/* process error vector */

	KalMat		P_k;	/* Process covariance matrix */
	KalMat		Q_k;	/* Process cov. matrix prediction error */
	KalMat		A_k;	/* State prediction matrix */
	KalMat		B_k;	/* Control matrix */
};

/*
 * A bunch of utility copy-in and copy-out functions to convert between the
 * C-like kalman_mat_t and kalman_vec_t types and the Eigen-library types.
 */
static void
kalmat_copyin(KalMat &mat, const kalman_mat_t *c_mat, unsigned dim)
{
	for (unsigned col = 0; col < dim; col++) {
		for (unsigned row = 0; row < dim; row++)
			mat(row, col) = KALMAN_MATxy(*c_mat, col, row);
	}
}

static void
kalmat_copyout(kalman_mat_t *c_mat, const KalMat &mat, unsigned dim)
{
	memset(c_mat, 0, sizeof (*c_mat));
	for (unsigned col = 0; col < dim; col++) {
		for (unsigned row = 0; row < dim; row++)
			KALMAN_MATxy(*c_mat, row, col) = mat(row, col);
	}
}

static void
kalvec_copyin(KalVec &vec, const kalman_vec_t *c_vec, unsigned dim)
{
	for (unsigned row = 0; row < dim; row++)
		vec(row) = c_vec->v[row];
}

static void
kalvec_copyout(kalman_vec_t *c_vec, const KalVec &vec, unsigned dim)
{
	memset(c_vec, 0, sizeof (*c_vec));
	for (unsigned row = 0; row < dim; row++)
		c_vec->v[row] = vec(row);
}

static void
kalmat_init_value(KalMat &mat, unsigned dim, kalman_real_t value)
{
	for (unsigned row = 0; row < dim; row++) {
		for (unsigned col = 0; col < dim; col++)
			mat(row, col) = value;
	}
}

static void
kalvec_init_value(KalVec &vec, unsigned dim, kalman_real_t value)
{
	for (unsigned row = 0; row < dim; row++)
		vec(row) = value;
}

/*
 * Allocates and returns a new Kalman filter. Use `kalman_free' to deallocate
 * the filter and all its associated resources.
 *
 * @param state_len Length of the filter's state vector (i.e. how many
 *	variables you want the filter to work on). This must be greater
 *	than zero and equal to or less than KALMAN_VEC_LEN (9 by default).
 *	If you need a Kalman filter that can track more variables, change
 *	the KALMAN_VEC_LEN macro defined in `kalman.h'.
 */
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
	kal->x_k = KalVec(state_len);
	kalvec_init_value(kal->x_k, state_len, NAN);
	kal->u_k = KalVec::Zero(state_len);
	kal->w_k = KalVec::Zero(state_len);

	kal->P_k = KalMat(state_len, state_len);
	kalmat_init_value(kal->P_k, state_len, NAN);
	kal->Q_k = KalMat::Zero(state_len, state_len);
	kal->A_k = KalMat(state_len, state_len);
	kalmat_init_value(kal->A_k, state_len, NAN);
	kal->B_k = KalMat::Zero(state_len, state_len);

	return (kal);
}

/*
 * Destroys a Kalman filter previous allocated using kalman_alloc.
 */
void
kalman_free(kalman_t *kal)
{
	delete kal;
}

/*
 * Returns the state vector length of an allocated Kalman filter. This is
 * also equal to the number of rows and columns in the filter's matrices.
 */
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
	kalvec_copyin(kal->x_k, state, kal->state_len);
}

/*
 * Returns the Kalman filter's current state vector.
 */
kalman_vec_t
kalman_get_state(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	kalman_vec_t v;
	kalvec_copyout(&v, kal->x_k, kal->state_len);
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
	kalvec_copyin(kal->u_k, control, kal->state_len);
}

/*
 * Returns the Kalman filter's current control vector.
 */
kalman_vec_t
kalman_get_cont(const kalman_t *kal)
{
	kalman_vec_t v;
	KAL_ASSERT(kal != NULL);
	kalvec_copyout(&v, kal->u_k, kal->state_len);
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
kalman_set_proc_noise(kalman_t *kal, const kalman_vec_t *proc_err)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(proc_err != NULL);
	kalvec_copyin(kal->w_k, proc_err, kal->state_len);
}

/*
 * Returns the Kalman filter's process error vector.
 */
kalman_vec_t
kalman_get_proc_noise(const kalman_t *kal)
{
	kalman_vec_t v;
	KAL_ASSERT(kal != NULL);
	kalvec_copyout(&v, kal->w_k, kal->state_len);
	return (v);
}

/*
 * Sets the Kalman filter's current covariance matrix (P).
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
	kalmat_copyin(kal->P_k, cov_mat, kal->state_len);
}

/*
 * Returns the Kalman filter's current covariance matrix.
 */
kalman_mat_t
kalman_get_cov_mat(const kalman_t *kal)
{
	kalman_mat_t m;
	KAL_ASSERT(kal != NULL);
	kalmat_copyout(&m, kal->P_k, kal->state_len);
	return (m);
}

/*
 * Sets the Kalman filter's process noise covariance matrix (Q).
 */
void
kalman_set_proc_noise_cov(kalman_t *kal, const kalman_mat_t *cov_mat_err)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(cov_mat_err != NULL);
	kalmat_copyin(kal->Q_k, cov_mat_err, kal->state_len);
}

/*
 * Returns the Kalman filter's process noise covariance matrix.
 */
kalman_mat_t
kalman_get_proc_noise_cov(const kalman_t *kal)
{
	kalman_mat_t m;
	KAL_ASSERT(kal != NULL);
	kalmat_copyout(&m, kal->Q_k, kal->state_len);
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
	kalmat_copyin(kal->A_k, pred_mat, kal->state_len);
}

/*
 * Returns the Kalman filter's current prediction matrix.
 */
kalman_mat_t
kalman_get_pred_mat(const kalman_t *kal)
{
	kalman_mat_t m;
	KAL_ASSERT(kal != NULL);
	kalmat_copyout(&m, kal->A_k, kal->state_len);
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
	kalmat_copyin(kal->B_k, cont_mat, kal->state_len);
}

/*
 * Returns the Kalman filter's control matrix.
 */
kalman_mat_t
kalman_get_cont_mat(const kalman_t *kal)
{
	kalman_mat_t m;
	KAL_ASSERT(kal != NULL);
	kalmat_copyout(&m, kal->B_k, kal->state_len);
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

	KalVec m(kal->state_len);
	KalMat m_cov_mat(kal->state_len, kal->state_len);
	KalMat obsv_model(kal->state_len, kal->state_len);
	KalMat obsv_model_T(kal->state_len, kal->state_len);
	KalVec x_k_pred(kal->state_len);
	KalMat P_k_pred(kal->state_len, kal->state_len);
	KalMat tmpmat(kal->state_len, kal->state_len);
	KalMat K(kal->state_len, kal->state_len);

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
	CHECK_MAT(kal->A_k);
	CHECK_MAT(kal->x_k);
	CHECK_MAT(kal->B_k);
	CHECK_MAT(kal->u_k);
	CHECK_MAT(kal->w_k);
	x_k_pred = (kal->A_k * kal->x_k) + (kal->B_k * kal->u_k) + kal->w_k;
	P_k_pred = (kal->A_k * kal->P_k * kal->A_k.transpose()) + kal->Q_k;
	CHECK_MAT(x_k_pred);
	CHECK_MAT(P_k_pred);

	if (measurement == NULL) {
		KAL_ASSERT(measurement_cov_mat == NULL);
		kal->x_k = x_k_pred;
		kal->P_k = P_k_pred;
		return;
	}
	KAL_ASSERT(!KALMAN_IS_NULL_VEC(*measurement));
	KAL_ASSERT(!KALMAN_IS_NULL_MAT(*measurement_cov_mat));
	KAL_ASSERT(measurement_cov_mat != NULL);

	kalvec_copyin(m, measurement, kal->state_len);
	CHECK_MAT(m);
	kalmat_copyin(m_cov_mat, measurement_cov_mat, kal->state_len);
	CHECK_MAT(m_cov_mat);

	if (observation_model_p != NULL) {
		kalmat_copyin(obsv_model, observation_model_p, kal->state_len);
		CHECK_MAT(obsv_model);
	} else {
		obsv_model.setIdentity();
	}
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
	CHECK_MAT(K);
	/*
	 * Update the estimate via the measurement:
	 *                /            \
	 * x  = x' + K * ( z  - H * x'  )
	 *  k    k        \ k    k   k /
	 */
	kal->x_k = x_k_pred + K * (m - obsv_model * x_k_pred);
	CHECK_MAT(kal->x_k);
	/*
	 * P  = P' - K * H  * P'
	 *  k    k        k    k
	 */
	kal->P_k = P_k_pred - K * obsv_model * P_k_pred;
	CHECK_MAT(kal->P_k);

	KAL_ASSERT3F(kal->P_k(0, 0), >=, 0);
	KAL_ASSERT(isfinite(kal->x_k(0)));
}

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

	KalVec m0(state_len), m1(state_len), m(state_len);
	KalMat cov0(state_len, state_len), cov1(state_len, state_len);
	KalMat cov(state_len, state_len);
	KalMat K(state_len, state_len);

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

	kalvec_copyin(m0, m0_in, state_len);
	kalmat_copyin(cov0, cov0_in, state_len);
	kalvec_copyin(m1, m1_in, state_len);
	kalmat_copyin(cov1, cov1_in, state_len);

	K = cov0 * (cov0 + cov1).inverse();
	/*
	 * m' = m0 + K(m1 - m0)
	 */
	m = m0 + K * (m1 - m0);
	kalvec_copyout(m_out, m, state_len);
	KAL_ASSERT(!KALMAN_IS_NULL_VEC(*m_out));
	/*
	 * var' = var0 - K * var0
	 */
	cov = cov0 - K * cov0;
	kalmat_copyout(cov_out, cov, state_len);
	KAL_ASSERT(!KALMAN_IS_NULL_MAT(*cov_out));
}

void
kalman_print_vec(const char *name, const kalman_vec_t *vec, unsigned state_len)
{
	printf("%s = (", name);
	for (unsigned i = 0; i < state_len; i++) {
#ifdef	KALMAN_REAL_LONG_DOUBLE
		printf("%11.4Lf%s", KALMAN_VECi(*vec, i),
		    (i + 1 < state_len) ? " " : "");
#else	/* !defined(KALMAN_REAL_LONG_DOUBLE) */
		printf("%11.4f%s", KALMAN_VECi(*vec, i),
		    (i + 1 < state_len) ? " " : "");
#endif	/* !defined(KALMAN_REAL_LONG_DOUBLE) */
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
#ifdef	KALMAN_REAL_LONG_DOUBLE
			printf("%11.4Lf%s", KALMAN_MATxy(*mat, c, r),
			    (c + 1 < state_len) ? " " : "");
#else	/* !defined(KALMAN_REAL_LONG_DOUBLE) */
			printf("%11.4f%s", KALMAN_MATxy(*mat, c, r),
			    (c + 1 < state_len) ? " " : "");
#endif	/* !defined(KALMAN_REAL_LONG_DOUBLE) */
		}
		printf(")\n");
	}
}

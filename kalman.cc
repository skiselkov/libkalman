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
#include "kalman.hh"

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

	kalman::matrix	x_k;	/* state vector */
	kalman::matrix	u_k;	/* control vector */
	kalman::matrix	w_k;	/* process error vector */

	kalman::matrix	P_k;	/* Process covariance matrix */
	kalman::matrix	Q_k;	/* Process cov. matrix prediction error */
	kalman::matrix	A_k;	/* State prediction matrix */
	kalman::matrix	B_k;	/* Control matrix */
};

/*
 * A bunch of utility copy-in and copy-out functions to convert between the
 * C-like kalman_(d)mat_t and kalman_vec_t types and the Eigen-library types.
 */
static void
mat_copyin(kalman::matrix &mat, const kalman_mat_t *c_mat, unsigned dim)
{
	for (unsigned col = 0; col < dim; col++) {
		for (unsigned row = 0; row < dim; row++)
			mat(row, col) = KALMAN_MATyx(*c_mat, row, col);
	}
}

static void
mat_copyout(kalman_mat_t *c_mat, const kalman::matrix &mat, unsigned dim)
{
	memset(c_mat, 0, sizeof (*c_mat));
	for (unsigned col = 0; col < dim; col++) {
		for (unsigned row = 0; row < dim; row++)
			KALMAN_MATyx(*c_mat, row, col) = mat(row, col);
	}
}

static void
dmat_copyin(kalman::matrix &mat, const kalman_dmat_t *dmat)
{
	KAL_ASSERT(dmat != NULL);
	KAL_ASSERT(dmat->cols == mat.cols());
	KAL_ASSERT(dmat->rows == mat.rows());
	for (unsigned row = 0; row < dmat->rows; row++) {
		for (unsigned col = 0; col < dmat->cols; col++)
			mat(row, col) = KAL_DMATyx(*dmat, row, col);
	}
}

static kalman_dmat_t *
dmat_copyout(const kalman::matrix &mat)
{
	kalman_dmat_t *dmat = kalman_dmat_alloc(mat.rows(), mat.cols());
	for (unsigned row = 0; row < dmat->rows; row++) {
		for (unsigned col = 0; col < dmat->cols; col++)
			KAL_DMATyx(*dmat, row, col) = mat(row, col);
	}
	return (dmat);
}

static void
vec_copyin(kalman::matrix &vec, const kalman_vec_t *c_vec, unsigned dim)
{
	KAL_ASSERT(vec.cols() == 1);
	for (unsigned row = 0; row < dim; row++)
		vec(row) = c_vec->v[row];
}

static void
vec_copyout(kalman_vec_t *c_vec, const kalman::matrix &vec, unsigned dim)
{
	memset(c_vec, 0, sizeof (*c_vec));
	KAL_ASSERT(vec.cols() == 1);
	for (unsigned row = 0; row < dim; row++)
		c_vec->v[row] = vec(row);
}

static void
mat_init_value(kalman::matrix &mat, kalman_real_t value)
{
	for (unsigned row = 0; row < mat.rows(); row++) {
		for (unsigned col = 0; col < mat.cols(); col++)
			mat(row, col) = value;
	}
}

/**
 * Allocates and returns a new kalman filter. Use `kalman_free' to deallocate
 * the filter and all its associated resources.
 *
 * @param state_len Length of the filter's state vector (i.e. how many
 *	variables you want the filter to work on). This must be greater
 *	than zero and equal to or less than KALMAN_VEC_LEN (9 by default).
 *	If you need a kalman filter that can track more variables, change
 *	the KALMAN_VEC_LEN macro defined in `kalman.h'.
 */
kalman_t *
kalman_alloc(unsigned state_len)
{
	kalman_t *kal = new kalman_t;

	KAL_ASSERT(state_len != 0);
	memset(kal, 0, sizeof (*kal));
	kal->state_len = state_len;
	/*
	 * We init the state vector, covariance matrix and prediction matrix
	 * to NAN values. This forces the caller to set them before calling
	 * kalman_step for the first time.
	 */
	kal->x_k = kalman::matrix(state_len, 1);
	mat_init_value(kal->x_k, NAN);
	kal->u_k = kalman::matrix::Zero(state_len, 1);
	kal->w_k = kalman::matrix::Zero(state_len, 1);

	kal->P_k = kalman::matrix(state_len, state_len);
	mat_init_value(kal->P_k, NAN);
	kal->Q_k = kalman::matrix::Zero(state_len, state_len);
	kal->A_k = kalman::matrix(state_len, state_len);
	mat_init_value(kal->A_k, NAN);
	kal->B_k = kalman::matrix::Zero(state_len, state_len);

	return (kal);
}

/**
 * Destroys a kalman filter previous allocated using kalman_alloc.
 */
void
kalman_free(kalman_t *kal)
{
	delete kal;
}

/**
 * Returns the state vector length of an allocated kalman filter. This is
 * also equal to the number of rows and columns in the filter's matrices.
 */
unsigned
kalman_get_state_len(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (kal->state_len);
}

/**
 * Sets the current state vector of the kalman filter.
 * Use this to set an initial state.
 *
 * Theoretical background:
 * The state vector x_k represents the state of the system we want to model
 * in the kalman filter. Each element of the vector represents one state
 * parameter. The state's parameters are updated according to the measurements
 * taken and the relationship between their covariance matrices.
 */
void
kalman_set_state(kalman_t *kal, const kalman_vec_t *state)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(state != NULL);
	vec_copyin(kal->x_k, state, kal->state_len);
}

/**
 * Returns the filter's current state vector.
 */
kalman_vec_t
kalman_get_state(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	kalman_vec_t v;
	vec_copyout(&v, kal->x_k, kal->state_len);
	return (v);
}

/**
 * Same as kalman_set_state, but uses a dynamically sizeable matrix.
 */
void
kalman_set_dstate(kalman_t *kal, const kalman_dmat_t *state)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(state != NULL);
	dmat_copyin(kal->x_k, state);
}

/**
 * Same as kalman_get_state, but returns a dynamically sizeable matrix.
 * The caller is responsible for freeing the returned matrix using free().
 */
kalman_dmat_t *
kalman_get_dstate(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (dmat_copyout(kal->x_k));
}

/**
 * Sets the current state vector of the kalman filter.
 * Use this to set an initial state.
 *
 * Theoretical background:
 * The state vector x_k represents the state of the system we want to model
 * in the kalman filter. Each element of the vector represents one state
 * parameter. The state's parameters are updated according to the measurements
 * taken and the relationship between their covariance matrices.
 */
void
kalman::set_state(kalman_t *kal, const kalman::matrix &state)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(state.rows() == kal->state_len);
	KAL_ASSERT(state.cols() == 1);
	kal->x_k = state;
}

/**
 * Returns the filter's current state vector.
 */
kalman::matrix
kalman::get_state(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (kal->x_k);
}

/**
 * Sets the filter's control vector.
 *
 * Theoretical background:
 * The control vector u_k is used to modify the prediction for the next
 * state. For example, if the kalman filter is used to track position
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
	vec_copyin(kal->u_k, control, kal->state_len);
}

/**
 * Returns the filter's current control vector.
 */
kalman_vec_t
kalman_get_cont(const kalman_t *kal)
{
	kalman_vec_t v;
	KAL_ASSERT(kal != NULL);
	vec_copyout(&v, kal->u_k, kal->state_len);
	return (v);
}

/**
 * Same as kalman_set_cont, but uses a dynamically sizeable matrix.
 */
void
kalman_set_dcont(kalman_t *kal, const kalman_dmat_t *control)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(control != NULL);
	dmat_copyin(kal->u_k, control);
}

/**
 * Same as kalman_get_cont, but returns a dynamically sizeable matrix.
 * The caller is responsible for freeing the returned matrix using free().
 */
kalman_dmat_t *
kalman_get_dcont(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (dmat_copyout(kal->u_k));
}

/**
 * Sets the filter's control vector.
 *
 * Theoretical background:
 * The control vector u_k is used to modify the prediction for the next
 * state. For example, if the kalman filter is used to track position
 * and velocity, and we have additional knowledge of the instantaneous
 * acceleration, we can supply the acceleration component in the control
 * vector and use it to modify the predicted position and velocity at the
 * next step.
 */
void
kalman::set_cont(kalman_t *kal, const kalman::matrix &cont)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(cont.rows() == kal->state_len);
	KAL_ASSERT(cont.cols() == 1);
	kal->u_k = cont;
}

/**
 * Returns the filter's current control vector.
 */
kalman::matrix
kalman::get_cont(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (kal->u_k);
}

/**
 * Sets the filter's process noise vector.
 *
 * Theoretical background:
 * The process error vector w_k is a value added on top of the predicted
 * state vector to generate the predicted next state vector. It can
 * represent something like a constant error accumulated in the process
 * of predicting the next state of the system.
 */
void
kalman_set_proc_noise(kalman_t *kal, const kalman_vec_t *proc_noise)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(proc_noise != NULL);
	vec_copyin(kal->w_k, proc_noise, kal->state_len);
}

/**
 * Returns the filter's process noise vector.
 */
kalman_vec_t
kalman_get_proc_noise(const kalman_t *kal)
{
	kalman_vec_t v;
	KAL_ASSERT(kal != NULL);
	vec_copyout(&v, kal->w_k, kal->state_len);
	return (v);
}

/**
 * Same as kalman_set_proc_noise, but uses a dynamically sizeable matrix.
 */
void
kalman_set_dproc_noise(kalman_t *kal, const kalman_dmat_t *proc_err)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(proc_err != NULL);
	dmat_copyin(kal->w_k, proc_err);
}

/**
 * Same as kalman_get_proc_noise, but returns a dynamically sizeable matrix.
 * The caller is responsible for freeing the returned matrix using free().
 */
kalman_dmat_t *
kalman_get_dproc_noise(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (dmat_copyout(kal->w_k));
}

/**
 * Sets the filter's process error vector.
 *
 * Theoretical background:
 * The process error vector w_k is a value added on top of the predicted
 * state vector to generate the predicted next state vector. It can
 * represent something like a constant error accumulated in the process
 * of predicting the next state of the system.
 */
void
kalman::set_proc_noise(kalman_t *kal, const kalman::matrix &proc_noise)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(proc_noise.rows() == kal->state_len);
	KAL_ASSERT(proc_noise.cols() == 1);
	kal->w_k = proc_noise;
}

/**
 * Returns the filter's process noise vector.
 */
kalman::matrix
kalman::get_proc_noise(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (kal->w_k);
}

/**
 * Sets the filter's current covariance matrix (P).
 *
 * Theoretical background:
 * The process covariance matrix P_k represents the uncertainty
 * (variance-covariance) in the current state of the system. You should
 * set an initial process covariance matrix to give the kalman filter a
 * reasonable starting point. The filter will then update the covariance
 * matrix every step, depending on the covariance matrices of the
 * prediction step and measurement.
 */
void
kalman_set_cov_mat(kalman_t *kal, const kalman_mat_t *cov_mat)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(cov_mat != NULL);
	mat_copyin(kal->P_k, cov_mat, kal->state_len);
}

/**
 * Returns the filter's current covariance matrix.
 */
kalman_mat_t
kalman_get_cov_mat(const kalman_t *kal)
{
	kalman_mat_t m;
	KAL_ASSERT(kal != NULL);
	mat_copyout(&m, kal->P_k, kal->state_len);
	return (m);
}

/*
 * Same as kalman_set_cov_mat, but uses a dynamically sizeable matrix.
 */
void
kalman_set_dcov_mat(kalman_t *kal, const kalman_dmat_t *cov_mat)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(cov_mat != NULL);
	dmat_copyin(kal->P_k, cov_mat);
}

/*
 * Same as kalman_get_cov_mat, but returns a dynamically sizeable matrix.
 * The caller is responsible for freeing the returned matrix using free().
 */
kalman_dmat_t *
kalman_get_dcov_mat(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (dmat_copyout(kal->P_k));
}

/**
 * Sets the filter's current covariance matrix (P).
 *
 * Theoretical background:
 * The process covariance matrix P_k represents the uncertainty
 * (variance-covariance) in the current state of the system. You should
 * set an initial process covariance matrix to give the kalman filter a
 * reasonable starting point. The filter will then update the covariance
 * matrix every step, depending on the covariance matrices of the
 * prediction step and measurement.
 */
void
kalman::set_cov_mat(kalman_t *kal, const kalman::matrix &cov_mat)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(cov_mat.rows() == kal->state_len);
	KAL_ASSERT(cov_mat.cols() == kal->state_len);
	kal->P_k = cov_mat;
}

/**
 * Returns the filter's current covariance matrix.
 */
kalman::matrix
kalman::get_cov_mat(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (kal->P_k);
}

/**
 * Sets the filter's process noise covariance matrix (Q).
 */
void
kalman_set_proc_noise_cov(kalman_t *kal, const kalman_mat_t *proc_noise_cov)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(proc_noise_cov != NULL);
	mat_copyin(kal->Q_k, proc_noise_cov, kal->state_len);
}

/**
 * Returns the filter's process noise covariance matrix.
 */
kalman_mat_t
kalman_get_proc_noise_cov(const kalman_t *kal)
{
	kalman_mat_t m;
	KAL_ASSERT(kal != NULL);
	mat_copyout(&m, kal->Q_k, kal->state_len);
	return (m);
}

/*
 * Same as kalman_set_proc_noise_cov, but uses a dynamically sizeable matrix.
 */
void
kalman_set_dproc_noise_cov(kalman_t *kal, const kalman_dmat_t *proc_noise_cov)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(proc_noise_cov != NULL);
	dmat_copyin(kal->Q_k, proc_noise_cov);
}

/*
 * Same as kalman_get_proc_noise_cov, but returns a dynamically sizeable matrix.
 * The caller is responsible for freeing the returned matrix using free().
 */
kalman_dmat_t *
kalman_get_dproc_noise_cov(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (dmat_copyout(kal->Q_k));
}

/**
 * Sets the filter's process noise covariance matrix (Q).
 */
void
kalman::set_proc_noise_cov(kalman_t *kal, const kalman::matrix &proc_noise_cov)
{
	KAL_ASSERT(kal != NULL);
	kal->Q_k = proc_noise_cov;
}

/**
 * Returns the filter's process noise covariance matrix.
 */
kalman::matrix
kalman::get_proc_noise_cov(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (kal->Q_k);
}

/*
 * Sets the filter's prediction matrix.
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
	mat_copyin(kal->A_k, pred_mat, kal->state_len);
}

/*
 * Returns the filter's current prediction matrix.
 */
kalman_mat_t
kalman_get_pred_mat(const kalman_t *kal)
{
	kalman_mat_t m;
	KAL_ASSERT(kal != NULL);
	mat_copyout(&m, kal->A_k, kal->state_len);
	return (m);
}

/*
 * Same as kalman_set_pred_mat, but uses a dynamically sizeable matrix.
 */
void
kalman_set_dpred_mat(kalman_t *kal, const kalman_dmat_t *pred_mat)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(pred_mat != NULL);
	dmat_copyin(kal->A_k, pred_mat);
}

/*
 * Same as kalman_get_pred_mat, but returns a dynamically sizeable matrix.
 * The caller is responsible for freeing the returned matrix using free().
 */
kalman_dmat_t *
kalman_get_dpred_mat(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (dmat_copyout(kal->A_k));
}

/*
 * Sets the filter's prediction matrix.
 *
 * Theoretical background:
 * This matrix is used to evolve the state of the system to the next
 * step. It is used to multiply the current state vector to arrive at
 * the next predicted system state.
 */
void
kalman::set_pred_mat(kalman_t *kal, const kalman::matrix &pred_mat)
{
	KAL_ASSERT(kal != NULL);
	kal->A_k = pred_mat;
}

/*
 * Returns the filter's current prediction matrix.
 */
kalman::matrix
kalman::get_pred_mat(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (kal->A_k);
}

/**
 * Sets the filter's control matrix.
 *
 * Theoretical background:
 * The control matrix adapts the control vector and transforms the
 * values in the control vector to apply properly to the corresponding
 * state vector components. For example, in the case of a position-velocity
 * kalman filter with acceleration in the control vector, the control matrix
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
	mat_copyin(kal->B_k, cont_mat, kal->state_len);
}

/**
 * Returns the filter's control matrix.
 */
kalman_mat_t
kalman_get_cont_mat(const kalman_t *kal)
{
	kalman_mat_t m;
	KAL_ASSERT(kal != NULL);
	mat_copyout(&m, kal->B_k, kal->state_len);
	return (m);
}

/**
 * Same as kalman_set_cont_mat, but uses a dynamically sizeable matrix.
 */
void
kalman_set_dcont_mat(kalman_t *kal, const kalman_dmat_t *cont_mat)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(cont_mat != NULL);
	dmat_copyin(kal->B_k, cont_mat);
}

/**
 * Same as kalman_get_cont_mat, but returns a dynamically sizeable matrix.
 * The caller is responsible for freeing the returned matrix using free().
 */
kalman_dmat_t *
kalman_get_dcont_mat(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (dmat_copyout(kal->B_k));
}

/**
 * Returns the filter's control matrix.
 */
void
kalman::set_cont_mat(kalman_t *kal, const kalman::matrix &cont_mat)
{
	KAL_ASSERT(kal != NULL);
	kal->B_k = cont_mat;
}

kalman::matrix
kalman::get_cont_mat(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
	return (kal->B_k);
}

/*
 * Performs one step of the kalman filter, updating the filter's internal
 * state vector and process covariance matrix with a new measurement and
 * measurement covariance matrix.
 *
 * Theoretical background:
 * The kalman filter combines the current state with a prediction and
 * measurement to arrive at the next state of the system. The state's
 * preference for either using its internal predicted state, or the
 * measured state depends on the ratio of the predicted process
 * covariance matrix (calculated internally by the filter) and the
 * measurement covariance matrix supplied in this call.
 */
void
kalman::step(kalman_t *kal, const kalman::matrix &m,
    const kalman::matrix &m_cov_mat, const kalman::matrix &obsv_model)
{
	KAL_ASSERT(kal != NULL);
	KAL_ASSERT(!isnan(kal->x_k(0)));
	KAL_ASSERT(!isnan(kal->P_k(0, 0)));
	KAL_ASSERT(!isnan(kal->A_k(0, 0)));

	kalman::matrix obsv_model_T(obsv_model.rows(), obsv_model.cols());
	kalman::matrix x_k_pred(kal->state_len, 1);
	kalman::matrix P_k_pred(kal->state_len, kal->state_len);
	kalman::matrix tmpmat(kal->state_len, kal->state_len);
	kalman::matrix K(kal->state_len, kal->state_len);

	obsv_model_T = obsv_model.transpose();

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

	CHECK_MAT(m);
	CHECK_MAT(m_cov_mat);
	/*
	 * Compute the kalman gain:
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

void
kalman::step(kalman_t *kal, const kalman::matrix &m,
    const kalman::matrix &m_cov_mat)
{
	KAL_ASSERT(kal != NULL);
	kalman::matrix obsv_model(kal->state_len, kal->state_len);
	obsv_model.setIdentity();
	kalman::step(kal, m, m_cov_mat, obsv_model);
}

/*
 * Performs a no-measurement update of the filter using only
 * the prediction logic. Use this function when you don't have
 * a new measurement to supply to the filter.
 */
void
kalman::step(kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);

	kal->x_k = (kal->A_k * kal->x_k) + (kal->B_k * kal->u_k) +
	    kal->w_k;
	kal->P_k = (kal->A_k * kal->P_k * kal->A_k.transpose()) + kal->Q_k;
}

/**
 * Performs one step of the kalman filter, updating the filter's internal
 * state vector and process covariance matrix with a new measurement and
 * measurement covariance matrix.
 *
 * @param measurement A vector containing the measurement to be integrated
 *	into the filter's model. You can pass NULL here - this performs a
 *	null step by just performing the prediction and establishes the
 *	prediction as the new filter state.
 * @param measurement_cov_mat The measurement covariance matrix. Must be
 *	present when a measurement is supplied, or NULL when no measurement
 *	is supplied.
 * @param observation_model An optional matrix (you may pass NULL) which is
 *	used to adapt the measurement vector to the shape of the filter's
 *	state vector and covariance matrix. Because the fixed-size matrix
 *	API only works with the matrices being square and the same dimension
 *	as the filter, this argument is pretty useless in this function. It
 *	is provided primarily for compatibility reasons. If you plan to use
 *	the observation model properly, use dynamically sized matrices and
 *	kalman_dstep instead.
 *
 * Theoretical background:
 * The kalman filter combines the current state with a prediction and
 * measurement to arrive at the next state of the system. The state's
 * preference for either using its internal predicted state, or the
 * measured state depends on the ratio of the predicted process
 * covariance matrix (calculated internally by the filter) and the
 * measurement covariance matrix supplied in this call.
 */
void
kalman_step(kalman_t *kal, const kalman_vec_t *measurement,
    const kalman_mat_t *measurement_cov_mat,
    const kalman_mat_t *observation_model)
{
	KAL_ASSERT(kal != NULL);

	if (measurement == NULL) {
		KAL_ASSERT(measurement_cov_mat == NULL);
		kalman::step(kal);
		return;
	}

	KAL_ASSERT(!KALMAN_IS_NULL_VEC(*measurement));
	KAL_ASSERT(!KALMAN_IS_NULL_MAT(*measurement_cov_mat));
	KAL_ASSERT(measurement_cov_mat != NULL);

	kalman::matrix m(kal->state_len, 1);
	kalman::matrix m_cov_mat(kal->state_len, kal->state_len);

	vec_copyin(m, measurement, kal->state_len);
	CHECK_MAT(m);
	mat_copyin(m_cov_mat, measurement_cov_mat, kal->state_len);
	CHECK_MAT(m_cov_mat);

	if (observation_model != NULL) {
		kalman::matrix obsv_model(kal->state_len, kal->state_len);

		mat_copyin(obsv_model, observation_model, kal->state_len);
		CHECK_MAT(obsv_model);
		kalman::step(kal, m, m_cov_mat, obsv_model);
	} else {
		kalman::step(kal, m, m_cov_mat);
	}
}

/**
 * Same as kalman_step, but uses dynamically sizeable matrices.
 */
void
kalman_dstep(kalman_t *kal, const kalman_dmat_t *measurement,
    const kalman_dmat_t *measurement_cov_mat,
    const kalman_dmat_t *observation_model)
{
	ASSERT(kal != NULL);

	if (measurement == NULL) {
		KAL_ASSERT(measurement_cov_mat == NULL);
		kalman::step(kal);
		return;
	}
	KAL_ASSERT(measurement_cov_mat != NULL);

	kalman::matrix m(measurement->rows, measurement->cols);
	kalman::matrix m_cov_mat(measurement_cov_mat->rows, measurement_cov_mat->cols);

	dmat_copyin(m, measurement);
	dmat_copyin(m_cov_mat, measurement_cov_mat);

	CHECK_MAT(m);
	CHECK_MAT(m_cov_mat);

	if (observation_model != NULL) {
		kalman::matrix obsv_model(observation_model->rows,
		    observation_model->cols);

		dmat_copyin(obsv_model, observation_model);
		CHECK_MAT(obsv_model);

		kalman::step(kal, m, m_cov_mat, obsv_model);
	} else {
		kalman::step(kal, m, m_cov_mat);
	}
}

/**
 * Takes two scalar measurements and combines them together according to
 * their variances (sigma^2), producing a new combined measurement and
 * variance. Use this if you want to update your kalman filter with
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
kalman_combine_s(kalman_real_t m0, kalman_real_t var0,
    kalman_real_t m1, kalman_real_t var1,
    kalman_real_t *m_out, kalman_real_t *var_out)
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

/**
 * Takes two measurements and their covariance matrices and combines them
 * into a single measurement according to their covariances. The resulting
 * measurement and covariance is suitable for passing in kalman::step. Use
 * this to implement sensor fusion where multiple independent sensors are
 * supplying measurement data that will update the filter in a single step.
 *
 * @param m0_in First measurement.
 * @param cov0_in Covariance of first measurement.
 * @param m1_in Second measurement.
 * @param cov1_in Covariance of second measurement.
 * @param m_out Output that will be populated with the combined measurement.
 * @param cov_out Output that will be populated with the combined covariance.
 */
void
kalman::combine(const kalman::matrix &m0, const kalman::matrix &cov0,
    const kalman::matrix &m1, const kalman::matrix &cov1,
    kalman::matrix &m_out, kalman::matrix &cov_out)
{
	KAL_ASSERT(m0.rows() == cov0.rows());
	KAL_ASSERT(m0.cols() == 1);
	KAL_ASSERT(cov0.rows() == cov0.cols());
	KAL_ASSERT(m1.rows() == cov1.rows());
	KAL_ASSERT(m1.cols() == 1);
	KAL_ASSERT(cov1.rows() == cov1.cols());
	KAL_ASSERT(cov0.rows() == cov1.rows());

	kalman::matrix K = cov0 * (cov0 + cov1).inverse();

	/*
	 * m' = m0 + K(m1 - m0)
	 */
	m_out = m0 + K * (m1 - m0);
	/*
	 * var' = var0 - K * var0
	 */
	cov_out = cov0 - K * cov0;
}

/**
 * Same as kalman_combine, but uses statically-sized vectors & matrices
 * for the measurements and covariances.
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
void
kalman_combine_v(unsigned state_len,
    const kalman_vec_t *m0_in, const kalman_mat_t *cov0_in,
    const kalman_vec_t *m1_in, const kalman_mat_t *cov1_in,
    kalman_vec_t *m_out, kalman_mat_t *cov_out)
{
	KAL_ASSERT(state_len != 0);
	KAL_ASSERT(state_len <= KALMAN_VEC_LEN);

	kalman::matrix m0(state_len, 1), m1(state_len, 1), m(state_len, 1);
	kalman::matrix cov0(state_len, state_len), cov1(state_len, state_len);
	kalman::matrix cov(state_len, state_len);
	kalman::matrix K(state_len, state_len);

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

	vec_copyin(m0, m0_in, state_len);
	mat_copyin(cov0, cov0_in, state_len);
	vec_copyin(m1, m1_in, state_len);
	mat_copyin(cov1, cov1_in, state_len);

	K = cov0 * (cov0 + cov1).inverse();
	/*
	 * m' = m0 + K(m1 - m0)
	 */
	m = m0 + K * (m1 - m0);
	vec_copyout(m_out, m, state_len);
	KAL_ASSERT(!KALMAN_IS_NULL_VEC(*m_out));
	/*
	 * var' = var0 - K * var0
	 */
	cov = cov0 - K * cov0;
	mat_copyout(cov_out, cov, state_len);
	KAL_ASSERT(!KALMAN_IS_NULL_MAT(*cov_out));
}

/**
 * Takes two measurements and their covariance matrices and combines them
 * into a single measurement according to their covariances. The resulting
 * measurement and covariance is suitable for passing in kalman::step. Use
 * this to implement sensor fusion where multiple independent sensors are
 * supplying measurement data that will update the filter in a single step.
 *
 * @param m0_in First measurement.
 * @param cov0_in Covariance of first measurement.
 * @param m1_in Second measurement.
 * @param cov1_in Covariance of second measurement.
 * @param m_out Output that will be populated with the combined measurement.
 *	The caller is responsible for freeing the returned matrix using free().
 * @param cov_out Output that will be populated with the combined covariance.
 *	The caller is responsible for freeing the returned matrix using free().
 */
void
kalman_combine_d(const kalman_dmat_t *m0_in, const kalman_dmat_t *cov0_in,
    const kalman_dmat_t *m1_in, const kalman_dmat_t *cov1_in,
    kalman_dmat_t **m_out, kalman_dmat_t **cov_out)
{
	KAL_ASSERT(m0_in != NULL);
	KAL_ASSERT(cov0_in != NULL);
	KAL_ASSERT(m1_in != NULL);
	KAL_ASSERT(cov1_in != NULL);
	KAL_ASSERT(m_out != NULL);
	KAL_ASSERT(cov_out != NULL);

	unsigned state_len = m0_in->rows;
	kalman::matrix m0(state_len, 1), m1(state_len, 1), m(state_len, 1);
	kalman::matrix cov0(state_len, state_len), cov1(state_len, state_len);
	kalman::matrix cov(state_len, state_len), K(state_len, state_len);

	dmat_copyin(m0, m0_in);
	dmat_copyin(cov0, cov0_in);
	dmat_copyin(m1, m1_in);
	dmat_copyin(cov1, cov1_in);

	K = cov0 * (cov0 + cov1).inverse();
	/*
	 * m' = m0 + K(m1 - m0)
	 */
	m = m0 + K * (m1 - m0);
	*m_out = dmat_copyout(m);
	/*
	 * var' = var0 - K * var0
	 */
	cov = cov0 - K * cov0;
	*cov_out = dmat_copyout(cov);
}

void
kalman_print_vec(const char *name, const kalman_vec_t *vec, unsigned state_len)
{
	KAL_ASSERT(name != NULL);
	KAL_ASSERT(vec != NULL);

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

	KAL_ASSERT(name != NULL);
	KAL_ASSERT(mat != NULL);

	printf("%s = (", name);
	memset(leadspace, ' ', len);
	leadspace[len] = '\0';
	leadspace[len - 1] = '(';
	for (unsigned r = 0; r < state_len; r++) {
		if (r > 0)
			printf("%s", leadspace);
		for (unsigned c = 0; c < state_len; c++) {
#ifdef	KALMAN_REAL_LONG_DOUBLE
			printf("%11.4Lf%s", KALMAN_MATyx(*mat, r, c),
			    (c + 1 < state_len) ? " " : "");
#else	/* !defined(KALMAN_REAL_LONG_DOUBLE) */
			printf("%11.4f%s", KALMAN_MATyx(*mat, r, c),
			    (c + 1 < state_len) ? " " : "");
#endif	/* !defined(KALMAN_REAL_LONG_DOUBLE) */
		}
		printf(")\n");
	}
}

void
kalman_print_dmat(const char *name, const kalman_dmat_t *mat)
{
	int len = strlen(name) + 4;
	char leadspace[len + 1];

	KAL_ASSERT(name != NULL);
	KAL_ASSERT(mat != NULL);

	printf("%s = (", name);
	memset(leadspace, ' ', len);
	leadspace[len] = '\0';
	leadspace[len - 1] = '(';
	for (unsigned r = 0; r < mat->rows; r++) {
		if (r > 0)
			printf("%s", leadspace);
		for (unsigned c = 0; c < mat->cols; c++) {
#ifdef	KALMAN_REAL_LONG_DOUBLE
			printf("%11.4Lf%s", KAL_DMATyx(*mat, r, c),
			    (c + 1 < mat->cols) ? " " : "");
#else	/* !defined(KALMAN_REAL_LONG_DOUBLE) */
			printf("%11.4f%s", KAL_DMATyx(*mat, r, c),
			    (c + 1 < mat->cols) ? " " : "");
#endif	/* !defined(KALMAN_REAL_LONG_DOUBLE) */
		}
		printf(")\n");
	}
}

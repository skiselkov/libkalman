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

#ifdef	KALMAN_USE_MATHC
#include "mathc.h"
CTASSERT(sizeof (mfloat_t) == sizeof (double));
#else	/* !defined(KALMAN_USE_MATHC) */
#include <matrix.h>
#include <matrix2.h>
CTASSERT(sizeof (Real) == sizeof (double));
#endif	/* !defined(KALMAN_USE_MATHC) */

#include "kalman.h"

#ifndef	KALMAN_USE_MATHC

#define	MAKE_MAT(kal)	m_get((kal)->state_len, (kal)->state_len)
#define	MAKE_VEC(kal)	v_get((kal)->state_len)

#define	KAL_MAT_COPYIN(mesch_mat, kal_mat) \
	do { \
		MAT *m_mat = (mesch_mat); \
		const kalman_mat_t *k_mat = (kal_mat); \
		for (unsigned col = 0; col < m_mat->n; col++) { \
			for (unsigned row = 0; row < m_mat->m; row++) { \
				m_set_val(m_mat, row, col, \
				    KALMAN_MATxy(*k_mat, col, row)); \
			} \
		} \
	} while (0)

#define	KAL_MAT_COPYOUT(kal_mat, mesch_mat) \
	do { \
		kalman_mat_t *k_mat = (kal_mat); \
		const MAT *m_mat = (mesch_mat); \
		for (unsigned col = 0; col < m_mat->n; col++) { \
			for (unsigned row = 0; row < m_mat->m; row++) { \
				KALMAN_MATxy(*k_mat, row, col) = \
				    m_get_val(m_mat, col, row); \
			} \
		} \
	} while (0)

#define	KAL_VEC_COPYIN(mesch_vec, kal_vec) \
	memcpy((mesch_vec)->ve, (kal_vec)->v, \
	    sizeof (double) * (mesch_vec)->dim)

#define	KAL_VEC_COPYOUT(kal_vec, mesch_vec) \
	memcpy((kal_vec)->v, (mesch_vec)->ve, \
	    sizeof (double) * (mesch_vec)->dim)

#endif	/* !defined(KALMAN_USE_MATHC) */

struct kalman_s {
	unsigned	state_len;

#ifdef	KALMAN_USE_MATHC
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
#else	/* !defined(KALMAN_USE_MATHC) */
	/* state vector */
	VEC		*x_k;
	/* control vector */
	VEC		*u_k;
	/* process error vector */
	VEC		*w_k;

	/* Process covariance matrix */
	MAT		*P_k;
	/* Process covariance matrix prediction error */
	MAT		*Q_k;
	/* State prediction matrix */
	MAT		*A_k;
	/* State prediction matrix - transposed */
	MAT		*A_k_T;
	/* Control matrix */
	MAT		*B_k;
	/* Sensor covariance matrix */
	MAT		*R_k;
	/*
	 * Runtime state - these get reused in kalman_step
	 */
	MAT		*tmpmat[3];
	MAT		*m_cov_mat;
	MAT		*observation_model;
	MAT		*observation_model_T;
	MAT		*P_k_pred;
	MAT		*K;

	VEC		*tmpvec[3];
	VEC		*x_k_pred;
	VEC		*m;
#endif	/* !defined(KALMAN_USE_MATHC) */
};

kalman_t *
kalman_alloc(unsigned state_len)
{
	kalman_t *kal = safe_calloc(1, sizeof (kalman_t));

	KAL_ASSERT(state_len != 0);
	KAL_ASSERT3U(state_len, <=, KALMAN_VEC_LEN);
	kal->state_len = state_len;

#ifdef	KALMAN_USE_MATHC
	/*
	 * We initialize the minimum set of vectors and matrices to null
	 * vector/matrix values, to force the user to set at least these
	 * parameters before invoking kalman_step for the first time.
	 */
	kal->x_k = KALMAN_NULL_VEC;
	kal->A_k = KALMAN_NULL_MAT;
	kal->P_k = KALMAN_NULL_MAT;
#else	/* !defined(KALMAN_USE_MATHC) */
	kal->x_k = MAKE_VEC(kal);
	kal->x_k->ve[0] = NAN;
	kal->u_k = MAKE_VEC(kal);
	kal->w_k = MAKE_VEC(kal);

	kal->P_k = MAKE_MAT(kal);
	kal->P_k->me[0][0] = NAN;
	kal->Q_k = MAKE_MAT(kal);
	kal->A_k = MAKE_MAT(kal);
	kal->A_k->me[0][0] = NAN;
	kal->A_k_T = MAKE_MAT(kal);
	kal->B_k = MAKE_MAT(kal);
	kal->R_k = MAKE_MAT(kal);

	for (int i = 0; i < 3; i++) {
		kal->tmpmat[i] = MAKE_MAT(kal);
		kal->tmpvec[i] = MAKE_VEC(kal);
	}

	kal->m_cov_mat = MAKE_MAT(kal);
	kal->observation_model = MAKE_MAT(kal);
	kal->observation_model_T = MAKE_MAT(kal);
	kal->P_k_pred = MAKE_MAT(kal);
	kal->K = MAKE_MAT(kal);

	kal->x_k_pred = MAKE_VEC(kal);
	kal->m = MAKE_VEC(kal);
#endif	/* !defined(KALMAN_USE_MATHC) */

	return (kal);
}

void
kalman_free(kalman_t *kal)
{
#ifndef	KALMAN_USE_MATHC
	V_FREE(kal->x_k);
	V_FREE(kal->u_k);
	V_FREE(kal->w_k);

	M_FREE(kal->P_k);
	M_FREE(kal->Q_k);
	M_FREE(kal->A_k);
	M_FREE(kal->A_k_T);
	M_FREE(kal->B_k);
	M_FREE(kal->R_k);

	for (int i = 0; i < 3; i++) {
		M_FREE(kal->tmpmat[i]);
		V_FREE(kal->tmpvec[i]);
	}

	M_FREE(kal->m_cov_mat);
	M_FREE(kal->observation_model);
	M_FREE(kal->observation_model_T);
	M_FREE(kal->P_k_pred);
	M_FREE(kal->K);

	V_FREE(kal->x_k_pred);
	V_FREE(kal->m);
#endif	/* !defined(KALMAN_USE_MATHC) */

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
#ifdef	KALMAN_USE_MATHC
	kal->x_k = *state;
#else
	KAL_VEC_COPYIN(kal->x_k, state);
#endif
}

/*
 * Returns the Kalman filter's current state vector.
 */
kalman_vec_t
kalman_get_state(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
#ifdef	KALMAN_USE_MATHC
	return (kal->x_k);
#else	/* !defined(KALMAN_USE_MATHC) */
	kalman_vec_t v;
	KAL_VEC_COPYOUT(&v, kal->x_k);
	return (v);
#endif	/* !defined(KALMAN_USE_MATHC) */
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
#ifdef	KALMAN_USE_MATHC
	kal->u_k = *control;
#else
	KAL_VEC_COPYIN(kal->u_k, control);
#endif
}

/*
 * Returns the Kalman filter's current control vector.
 */
kalman_vec_t
kalman_get_cont(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
#ifdef	KALMAN_USE_MATHC
	return (kal->u_k);
#else	/* !defined(KALMAN_USE_MATHC) */
	kalman_vec_t v;
	KAL_VEC_COPYOUT(&v, kal->u_k);
	return (v);
#endif	/* !defined(KALMAN_USE_MATHC) */
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
#ifdef	KALMAN_USE_MATHC
	kal->w_k = *proc_err;
#else
	KAL_VEC_COPYIN(kal->w_k, proc_err);
#endif
}

/*
 * Returns the Kalman filter's process error vector.
 */
kalman_vec_t
kalman_get_proc_err(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
#ifdef	KALMAN_USE_MATHC
	return (kal->w_k);
#else	/* !defined(KALMAN_USE_MATHC) */
	kalman_vec_t v;
	KAL_VEC_COPYOUT(&v, kal->w_k);
	return (v);
#endif	/* !defined(KALMAN_USE_MATHC) */
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
#ifdef	KALMAN_USE_MATHC
	kal->P_k = *cov_mat;
#else
	KAL_MAT_COPYIN(kal->P_k, cov_mat);
#endif
}

/*
 * Returns the Kalman filter's current process covariance matrix.
 */
kalman_mat_t
kalman_get_cov_mat(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
#ifdef	KALMAN_USE_MATHC
	return (kal->P_k);
#else	/* !defined(KALMAN_USE_MATHC) */
	kalman_mat_t m;
	KAL_MAT_COPYOUT(&m, kal->P_k);
	return (m);
#endif	/* !defined(KALMAN_USE_MATHC) */
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
#if	KALMAN_USE_MATHC
	kal->Q_k = *cov_mat_err;
#else
	KAL_MAT_COPYIN(kal->Q_k, cov_mat_err);
#endif
}

/*
 * Returns the Kalman filter's process covariance covariance matrix
 * prediction error.
 */
kalman_mat_t
kalman_get_cov_mat_err(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
#ifdef	KALMAN_USE_MATHC
	return (kal->Q_k);
#else	/* !defined(KALMAN_USE_MATHC) */
	kalman_mat_t m;
	KAL_MAT_COPYOUT(&m, kal->Q_k);
	return (m);
#endif	/* !defined(KALMAN_USE_MATHC) */
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
#ifdef	KALMAN_USE_MATHC
	kal->A_k = *pred_mat;
	mat4_transpose(kal->A_k_T.m, kal->A_k.m);
#else	/* !defined(KALMAN_USE_MATHC) */
	KAL_MAT_COPYIN(kal->A_k, pred_mat);
	m_transp(kal->A_k, kal->A_k_T);
#endif	/* !defined(KALMAN_USE_MATHC) */
}

/*
 * Returns the Kalman filter's current prediction matrix.
 */
kalman_mat_t
kalman_get_pred_mat(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
#ifdef	KALMAN_USE_MATHC
	return (kal->A_k);
#else	/* !defined(KALMAN_USE_MATHC) */
	kalman_mat_t m;
	KAL_MAT_COPYOUT(&m, kal->A_k);
	return (m);
#endif	/* !defined(KALMAN_USE_MATHC) */
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
#ifdef	KALMAN_USE_MATHC
	kal->B_k = *cont_mat;
#else
	KAL_MAT_COPYIN(kal->B_k, cont_mat);
#endif
}

/*
 * Returns the Kalman filter's control matrix.
 */
kalman_mat_t
kalman_get_cont_mat(const kalman_t *kal)
{
	KAL_ASSERT(kal != NULL);
#ifdef	KALMAN_USE_MATHC
	return (kal->B_k);
#else	/* !defined(KALMAN_USE_MATHC) */
	kalman_mat_t m;
	KAL_MAT_COPYOUT(&m, kal->B_k);
	return (m);
#endif	/* !defined(KALMAN_USE_MATHC) */
}

#ifdef	KALMAN_USE_MATHC

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

#endif	/* KALMAN_USE_MATHC */

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
	KAL_ASSERT(measurement != NULL);
	KAL_ASSERT(measurement_cov_mat != NULL);

#ifndef	KALMAN_USE_MATHC
	KAL_VEC_COPYIN(kal->m, measurement);
	KAL_MAT_COPYIN(kal->m_cov_mat, measurement_cov_mat);
	if (observation_model_p != NULL)
		KAL_MAT_COPYIN(kal->observation_model, observation_model_p);
	else
		m_ident(kal->observation_model);
	m_transp(kal->observation_model, kal->observation_model_T);

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
	mv_mlt(kal->A_k, kal->x_k, kal->tmpvec[0]);
	mv_mlt(kal->B_k, kal->u_k, kal->tmpvec[1]);
	v_add(kal->tmpvec[0], kal->tmpvec[1], kal->tmpvec[2]);
	v_add(kal->tmpvec[2], kal->w_k, kal->x_k_pred);

	m_mlt(kal->A_k, kal->P_k, kal->tmpmat[0]);
	m_mlt(kal->tmpmat[0], kal->A_k_T, kal->tmpmat[1]);
	m_add(kal->tmpmat[1], kal->Q_k, kal->P_k_pred);
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
	m_mlt(kal->P_k_pred, kal->observation_model_T, kal->tmpmat[0]);
	/* denominator portion */
	m_mlt(kal->observation_model, kal->tmpmat[0], kal->tmpmat[1]);
	m_add(kal->tmpmat[1], kal->m_cov_mat, kal->tmpmat[2]);
	m_inverse(kal->tmpmat[2], kal->tmpmat[1]);
	m_mlt(kal->tmpmat[0], kal->tmpmat[1], kal->K);
	/*
	 * Update the estimate via the measurement:
	 *                /            \
	 * x  = x' + K * ( z  - H * x'  )
	 *  k    k        \ k    k   k /
	 */
	mv_mlt(kal->observation_model, kal->x_k_pred, kal->tmpvec[0]);
	v_sub(kal->m, kal->tmpvec[0], kal->tmpvec[1]);
	mv_mlt(kal->K, kal->tmpvec[1], kal->tmpvec[0]);
	v_add(kal->x_k_pred, kal->tmpvec[0], kal->x_k);
	/*
	 * P  = P' - K * H  * P'
	 *  k    k        k    k
	 */
	m_mlt(kal->K, kal->observation_model, kal->tmpmat[0]);
	m_mlt(kal->tmpmat[0], kal->P_k_pred, kal->tmpmat[1]);
	m_sub(kal->P_k_pred, kal->tmpmat[1], kal->P_k);

	KAL_ASSERT3F(kal->P_k->me[0][0], >=, 0);
	KAL_ASSERT(isfinite(v_get_val(kal->x_k, 0)));
#else	/* KALMAN_USE_MATHC */
	kalman_vec_t x_k_pred, tmpvec1, tmpvec2;
	kalman_mat_t P_k_pred, tmpmat1, tmpmat2;
	kalman_mat_t K_numer, K_denom, K_denom_inv;
	kalman_mat_t K = KALMAN_ZERO_MAT;
	kalman_mat_t observation_model, observation_model_T;

	KAL_ASSERT(!KALMAN_IS_NULL_VEC(kal->x_k));
	KAL_ASSERT(!KALMAN_IS_NULL_MAT(kal->A_k));
	KAL_ASSERT(!KALMAN_IS_NULL_MAT(kal->P_k));

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
	vec4_multiply_mat4(tmpvec1.v, kal->x_k.v, kal->A_k.m);
	vec4_multiply_mat4(tmpvec2.v, kal->u_k.v, kal->B_k.m);
	vec4_add(x_k_pred.v, tmpvec1.v, tmpvec2.v);
	vec4_add(x_k_pred.v, x_k_pred.v, kal->w_k.v);

	mat4_multiply(tmpmat1.m, kal->A_k.m, kal->P_k.m);
	mat4_multiply(tmpmat2.m, tmpmat1.m, kal->A_k_T.m);
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
	mat4_multiply(tmpmat1.m, observation_model.m, K_numer.m);
	mat4_add(K_denom.m, tmpmat1.m, measurement_cov_mat->m);
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
	vec4_multiply_mat4(tmpvec1.v, x_k_pred.v, observation_model.m);
	vec4_subtract(tmpvec2.v, measurement->v, tmpvec1.v);
	vec4_multiply_mat4(tmpvec1.v, tmpvec2.v, K.m);
	vec4_add(kal->x_k.v, x_k_pred.v, tmpvec1.v);
	/*
	 * P  = P' - K * H  * P'
	 *  k    k        k    k
	 */
	mat4_multiply(tmpmat1.m, K.m, observation_model.m);
	mat4_multiply(tmpmat2.m, tmpmat1.m, P_k_pred.m);
	mat4_subtract(kal->P_k.m, P_k_pred.m, tmpmat2.m);

	KAL_ASSERT3F(kal->P_k.m[0], >=, 0);
	KAL_ASSERT(isfinite(KALMAN_VECi(kal->x_k, 0)));

#endif	/* KALMAN_USE_MATHC */
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
#ifndef	KALMAN_USE_MATHC
	/*
	 * x_k update
	 */
	mv_mlt(kal->A_k, kal->x_k, kal->tmpvec[0]);
	mv_mlt(kal->B_k, kal->u_k, kal->tmpvec[1]);
	v_add(kal->tmpvec[0], kal->tmpvec[1], kal->tmpvec[2]);
	v_add(kal->tmpvec[2], kal->w_k, kal->x_k);
	/*
	 * P_k update
	 */
	m_mlt(kal->A_k, kal->P_k, kal->tmpmat[0]);
	m_mlt(kal->tmpmat[0], kal->A_k_T, kal->tmpmat[1]);
	m_add(kal->tmpmat[1], kal->Q_k, kal->P_k);
#else	/* defined(KALMAN_USE_MATHC) */
	kalman_vec_t tmpvec1, tmpvec2, x_k;
	kalman_mat_t tmpmat1, tmpmat2;
	/*
	 * x_k update
	 */
	vec4_multiply_mat4(tmpvec1.v, kal->x_k.v, kal->A_k.m);
	vec4_multiply_mat4(tmpvec2.v, kal->u_k.v, kal->B_k.m);
	vec4_add(x_k.v, tmpvec1.v, tmpvec2.v);
	vec4_add(kal->x_k.v, x_k.v, kal->w_k.v);
	/*
	 * P_k update
	 */
	mat4_multiply(tmpmat1.m, kal->A_k.m, kal->P_k.m);
	mat4_multiply(tmpmat2.m, tmpmat1.m, kal->A_k_T.m);
	mat4_add(kal->P_k.m, tmpmat2.m, kal->Q_k.m);
#endif	/* defined(KALMAN_USE_MATHC) */
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
 * @param m0 First measurement.
 * @param var0 Covariance of first measurement.
 * @param m1 Second measurement.
 * @param var1 Covariance of second measurement.
 * @param m_out Output that will be populated with the combined measurement.
 * @param var_out Output that will be populated with the combined covariance.
 */
void kalman_combine_measurements_mat(unsigned state_len,
    const kalman_vec_t *restrict m0, const kalman_mat_t *restrict var0,
    const kalman_vec_t *restrict m1, const kalman_mat_t *restrict var1,
    kalman_vec_t *restrict m_out, kalman_mat_t *restrict var_out)
{
	VEC *v_m0, *v_m1, *v_X, *v_Y;
	MAT *m_var0, *m_var1, *m_X, *m_Y, *K;

	ASSERT(state_len != 0);
	ASSERT3U(state_len, <=, KALMAN_VEC_LEN);
	ASSERT(m0 != NULL);
	ASSERT(var0 != NULL);
	ASSERT(m1 != NULL);
	ASSERT(var1 != NULL);
	ASSERT(m_out != NULL);
	ASSERT(var_out != NULL);

	v_m0 = v_get(state_len);
	KAL_VEC_COPYIN(v_m0, m0);
	m_var0 = m_get(state_len, state_len);
	KAL_MAT_COPYIN(m_var0, var0);
	v_m1 = v_get(state_len);
	KAL_VEC_COPYIN(v_m1, m1);
	m_var1 = m_get(state_len, state_len);
	KAL_MAT_COPYIN(m_var1, var1);
	m_X = m_get(state_len, state_len);
	m_Y = m_get(state_len, state_len);
	K = m_get(state_len, state_len);
	v_X = v_get(state_len);
	v_Y = v_get(state_len);
	/*
	 * K = var0 / (var0 + var1)
	 */
	m_add(m_var0, m_var1, m_X);	/* X = var0 + var1 */
	m_inverse(m_X, m_Y);		/* Y = X^-1 */
	m_mlt(m_var0, m_Y, K);		/* K = var0 * Y */
	/*
	 * m' = m0 + K(m1 - m0)
	 */
	v_sub(v_m1, v_m0, v_X);		/* X = m1 - m0 */
	mv_mlt(K, v_X, v_Y);		/* Y = K * X */
	v_add(v_m0, v_Y, v_X);		/* X = m0 + Y */
	KAL_VEC_COPYOUT(m_out, v_X);
	/*
	 * var' = var0 - K * var0
	 */
	m_mlt(K, m_var0, m_X);		/* X = K * var0 */
	m_sub(m_var0, m_X, m_Y);	/* Y = var0 - X */
	KAL_MAT_COPYOUT(var_out, m_Y);

	V_FREE(v_m0);
	M_FREE(m_var0);
	V_FREE(v_m1);
	M_FREE(m_var1);

	M_FREE(m_X);
	M_FREE(m_Y);
	M_FREE(K);

	V_FREE(v_X);
	V_FREE(v_Y);
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

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

#ifndef	_LIBKALMAN_HH_
#define	_LIBKALMAN_HH_

#ifndef	__cplusplus
#error "Including this header requires using a C++ compiler"
#endif

#include "Eigen/Dense"
#include "kalman.h"

namespace kalman {

/**
 * The primary C++ interface type for the libkalman library. Use this type
 * to represent any vectors (rows=N, cols=1) and matrix (rows=N, cols=M)
 * data that you pass to the library via the C++ interface. The type is
 * fully RAII aware, so no special handling of memory allocation is
 * necessary.
 */
typedef Eigen::Matrix<kalman_real_t, Eigen::Dynamic, Eigen::Dynamic> matrix;

/*
 * Manipulating the filter's state vector and covariance matrix.
 */
void set_state(kalman_t *kal, const matrix &state);
matrix get_state(const kalman_t *kal);

void set_cov_mat(kalman_t *kal, const matrix &cov_mat);
matrix get_cov_mat(const kalman_t *kal);
/*
 * Manipulating the filter's state prediction matrix
 */
void set_pred_mat(kalman_t *kal, const matrix &pred_mat);
matrix get_pred_mat(const kalman_t *kal);
/*
 * Manipulating the filter's process noise vector and process noise covariance.
 */
void set_proc_noise(kalman_t *kal, const matrix &proc_noise);
matrix get_proc_noise(const kalman_t *kal);

void set_proc_noise_cov(kalman_t *kal, const matrix &proc_noise_cov);
matrix get_proc_noise_cov(const kalman_t *kal);
/*
 * Manipulating the filter's control vector and control matrix.
 */
void set_cont(kalman_t *kal, const matrix &cont);
matrix get_cont(const kalman_t *kal);

void set_cont_mat(kalman_t *kal, const matrix &cont_mat);
matrix get_cont_mat(const kalman_t *kal);
/*
 * Running the filter.
 */
void step(kalman_t *kal);
void step(kalman_t *kal, const matrix &m, const matrix &m_cov_mat);
void step(kalman_t *kal, const matrix &m, const matrix &m_cov_mat,
    const matrix &obsv_model);
/*
 * Utility functions.
 */
void combine(const matrix &m0, const matrix &cov0,
    const matrix &m1, const matrix &cov1, matrix &m_out, matrix &cov_out);

}

#endif	/* _LIBKALMAN_HH_ */

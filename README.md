# libkalman
Generic Kalman filter implementation with C and C++ interfaces.

## Usage

You simply need to add kalman.cc to your build process. The library
contains its own copy of the Eigen linear algebra library to do all the
mathematical heavy lifting. The library provides 3 different interfaces:

1) A simple (AND DEPRECATED!) C interface that uses statically-sized
   structures. These functions take fixed-size kalman_vec_t and
   kalman_mat_t arguments. If the filter needs fewer values than are in
   the structure, it will ignore the extraneous rows. Please note that
   the fixed-size structures are limited to 9-variable filters. If you
   need more variables in the filter, use the dynamic interface below.
   Please note that this interface is deprecated. Use the dynamically-
   sizeable C interface below for future development.
2) A C interface that uses dynamically-sized structures. These functions
   take kalman_dmat_t arguments. These are single-block allocation
   structures created using kalman_dmat_alloc. Once created using the
   ``kalman_dmat_alloc`` (and similar) functions, the matrix dimensions
   cannot be altered. Free the allocated matrix using the standard C
   library free() function.
   Dynamic interface functions have a 'd' in their name (e.g.
   ``kalman_set_dstate`` instead of ``kalman_set_state``).
3) A C++ interface that uses dynamically-sized matrices with RAII. These
   structures dynamically allocate memory as necessary and are aliases
   for the Eigen::Matrix type. All features of the Eigen::Matrix template
   type are available.
   The C++ interface functions are declared in ``kalman.hh`` and are
   located inside of the ``kalman`` namespace.

Should you need to, you can mix-and-match these interfaces. They all
operate on the same underlying state.

To construct a filter, use ``kalman_alloc`` with the number of state
variables you wish to track. Before being able to use the filter you
must initialize it:

* Setup the filter's initial state vector (x_k). Use the ``kalman_set_state``,
  ``kalman_set_dstate`` or ``kalman::set_state`` functions for this task.
* Setup the filter's initial covariance matrix (P_k). Use the
  ``kalman_set_cov_mat``, ``kalman_set_dcov_mat`` or ``kalman::set_cov_mat``
  functions for this task.
* Setup the filter's prediction matrix (A_k). Use the ``kalman_set_pred_mat``,
  ``kalman_set_dpred_mat`` or ``kalman::set_pred_mat`` functions for this task.

This is the minimum amount of setup required to get the filter working.
See ``kalman.h`` or ``kalman.hh`` for the other functions available to
initialize other useful properties of the filter (such as the control
vector, process noise, etc.). To evolve the filter by one step, use one
of the following functions:

```
kalman_step(kal, measurement, measurement_covariance, observation_model);
kalman_dstep(kal, measurement, measurement_covariance, observation_model);
kalman::step(kal, measurement, measurement_covariance, observation_model);
```

``measurement`` and ``measurement_covariance`` are the new measurement
vector and covariance matrix. If your measurement uses a different layout
than your filter's vector, you can pass the ``observation_model`` matrix,
which will be used to adapt the measurement to fit the filter's internal
state vector layout. If your measurement vector matches the filter's
state vector, simply pass NULL for the observation model (use the
2-argument version of kalman::step for the C++ interface).

In the C interfaces, you can also pass NULL for the ``measurement``
vector and ``measurement_covariance`` matrix. This will perform a
predictive-only (NULL) update of the filter. This is useful for cases
when your filter needs to be updated at fixed intervals, but your
measurement system doesn't always supply new measurements to update the
filter. Mathematically, this simply performs the predicted measurement
(x_k') and predicted covariance (P_k') calculations and then simply
installs these new predictions as the filter's current state vector and
covariance.

In the C++ interface, you can perform a NULL update by calling
``kalman::step`` with just the filter argument and no measurement.

You can retrieve the filter's current state vector and covariance matrix
using ``kalman_get_state`` and ``kalman_get_cov_mat``.

## Debugging

You will no doubt want to observe the internals of the filter during
development. The ``KALMAN_DEBUG_STATE`` and ``KALMAN_DEBUG_COV_MAT``
macros give you a quick shorthand that dumps the filter's current state
vector and covariance matrix to stdout. You can also use the more
generic ``KALMAN_DEBUG_VEC`` and ``KALMAN_DEBUG_MAT`` macros to extract
any other vector or matrix in the filter, e.g. to extract the control
vector:

```
kalman_t *kal;
...
KALMAN_DEBUG_VEC(kal, "cont_vect", kalman_get_cont);
```

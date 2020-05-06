# libkalman
Generic Kalman filter implementation with a pure C interface.

## Usage

You simply need to add kalman.cc to your build process. The library
contains its own copy of the Eigen linear algebra library to do all the
mathematical heavy lifting. libkalman provides its own C-compatible data
types for the matrices and vectors. Currently, the library is limited to
6th order Kalman filters (up to 6 tracked variables), but can be
trivially extended to support arbitrarily large states (see
``KALMAN_VEC_LEN`` in ``kalman.h``).

To construct a filter, use ``kalman_alloc`` with the number of state
variables you wish to track. Before being able to use the filter you need
to initialize its internal state:

*) Use ``kalman_set_state`` to initialize the filter's initial state vector
   (x_k).
*) Use ``kalman_set_cov_mat`` to initialize the filter's initial covariance
   matrix (P_k).
*) Use ``kalman_set_pred_mat`` to initialize the filter's prediction matrix
   (A_k).

This is the minimum amount of setup required to get the filter working.
See ``kalman.h`` for the other functions available to initialize other
useful properties of the filter (such as the control vector, process
error, etc.). You can then run the filter by invoking the ``kalman_step``
function:

```
kalman_step(kal, measurement, measurement_covariance, observation_model);
```

``measurement`` and ``measurement_covariance`` are the new measurement
vector and covariance matrix. If your measurement uses a different layout
than your filter's vector, you can pass the ``observation_model`` matrix,
which will be used to adapt the measurement to fit the filter's internal
state vector layout. If your measurement vector matches the filter's
state vector, simply pass NULL for the observation model.

You can also pass NULL for the ``measurement`` vector and
``measurement_covariance`` matrix. This will perform a predictive-only
(NULL) update of the filter. This is useful for cases when your filter
needs to be updated at fixed intervals, but your measurement system
doesn't always supply new measurements to update the filter.
Mathematically, this simply performs the predicted measurement (x_k') and
predicted covariance (P_k') calculations and then simply installs these
new predictions as the filter's current state vector and covariance.

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

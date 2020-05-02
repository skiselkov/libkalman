# libkalman
Generic Kalman filter implementation in C using the meschach or mathc
math library to do the heavy lifting. To use this library, simply
build it together with your code.

## Prerequisites

You have the option of either the Meschach library or MathC library
to provide the underlying matrix math engine.

1) Meschach: https://github.com/skiselkov/Meschach
   This is the preferred library to use with libkalman. It provides
   arbitrary size matrix support. For interfacing reasons, however,
   libkalman is limited to 6-dimensional Kalman filters, this can
   be trivially extended though, at the cost of more memory usage.

2) https://github.com/skiselkov/mathc/
   MathC is a very simple library that provides a bunch of useful
   mathematical operations besides matrices and it has a convenient
   API. The limitation for MathC is that it only supports to
   4x4 matrices, so libkalman will be limited to 4-dimensional filters.

To use Meschach, use the CMake file included in its repository. This
will give you a static library that can then be simply linked into your
resulting application. See XCompile.txt in the Meschach repository on
how to cross-compile for Windows from Linux.

To use MathC, include math.c and math.h in your build directly.
Also, be sure to define the following macros:

* MATHC_USE_DOUBLE_FLOATING_POINT - this enables double precision
  floating point support in MathC
* MATHC_USE_UNIONS - this requires at least rudimentary C99 support
  from your compiler and exposes MathC's data types as an anonymous
  union of an array accessor as well as explicit data fields.
* KALMAN_USE_MATHC - this configures libkalman to use MathC instead
  of Meschach (which is the default)

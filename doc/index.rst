.. jmtx documentation master file, created by
   sphinx-quickstart on Fri Jan 24 20:03:35 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

``jmtx`` documentation
======================

This is the documentation covering the :mod:`jmtx` library. The aim purpose of the
library is to provide a fast and clear implementation of iterative and direct solvers
for sparse matrices, ability to build and convert matrices to and from many sparse
formats, and perform algbraic operations with them.

The library is written in C99, which means it can be build by virtually any compiler.
All code is implemented for the four most commonly used C numeric real types:

- 32-bit floating point real numbers (``float``)
- 64-bit floating point real numbers (``double``)
- 32-bit floating point complex numbers (``_Complex float``)
- 64-bit floating point complex numbers (``_Complex double``)

For all of these, the library provides conversion functions between the complex and real
tipes of the same width, as well as the conversion between the two real and two complex
types.

.. mermaid::

   ---
   title: type conversion
   ---
   stateDiagram-v2

      float --> double : jmtxd_matrix_*_from_float
      double --> float : jmtx_matrix_*_from_double
      _Complex float --> _Complex double : jmtxz_matrix_*_from_cfloat
      _Complex double --> _Complex  float : jmtxc_matrix_*_from_cdouble
      float --> _Complex float : jmtxc_matrix_*_from_float
      _Complex float --> float : jmtx_matrix_*_from_cfloat
      double --> _Complex double :  jmtxz_matrix_*_from_double
      _Complex double --> double : jmtxd_matrix_*_from_cdouble

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   source/matrix_types

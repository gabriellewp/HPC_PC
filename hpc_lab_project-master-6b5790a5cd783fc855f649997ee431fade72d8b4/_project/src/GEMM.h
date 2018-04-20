#ifndef GEMM_H_
#define GEMM_H_

/// Generalized matrix-matrix multiplication for column-major storage
void DGEMM(unsigned M, unsigned N, unsigned K, double alpha, double const* A, unsigned ldA, double const* B, unsigned ldB, double beta, double* C, unsigned ldC);

/// Generalized matrix-matrix multiplication for column-major storage
/// A (from GlobalMatrices) transposed and padded
void opt_DGEMM1(unsigned M, unsigned N, unsigned K, double alpha, double const* A, unsigned ldA, double const* B, unsigned ldB, double beta, double* C, unsigned ldC);

/// Generalized matrix-matrix multiplication for column-major storage
/// for computeFlux() with padded tmp
void opt_DGEMM2(unsigned M, unsigned N, unsigned K, double alpha, double const* A, unsigned ldA, double const* B, unsigned ldB, double beta, double* C, unsigned ldC);

/// Generalized matrix-matrix multiplication for column-major storage
/// for sparse A as B
void opt_DGEMM4A(unsigned M, unsigned N, unsigned K, double alpha, double const* A, unsigned ldA, double const* B, unsigned ldB, double beta, double* C, unsigned ldC);

/// Generalized matrix-matrix multiplication for column-major storage
/// for sparse B as B
void opt_DGEMM4B(unsigned M, unsigned N, unsigned K, double alpha, double const* A, unsigned ldA, double const* B, unsigned ldB, double beta, double* C, unsigned ldC);

#endif // GEMM_H_

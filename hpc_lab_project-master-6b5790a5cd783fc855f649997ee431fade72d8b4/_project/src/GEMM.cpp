#include "GEMM.h"
#include "constants.h"
#include "omp.h"
#include <immintrin.h>

void DGEMM(unsigned M, unsigned N, unsigned K, double alpha, double const* A, unsigned ldA, double const* B, unsigned ldB, double beta, double* C, unsigned ldC) {
  for (unsigned j = 0; j < N; ++j) {
    for (unsigned i = 0; i < M; ++i) {
      double cij = 0.0;
      for (unsigned k = 0; k < K; ++k) {
        cij += A[k*ldA + i] * B[j*ldB + k];
      }
      C[j*ldC + i] = alpha * cij + beta * C[j*ldC + i];
    }
  }
}

#if !defined(PADDED_GLOBAL_MATRICES_LENGTH)

// same as DGEMM() for Fallback
void opt_DGEMM1(unsigned M, unsigned N, unsigned K, double alpha, double const* A, unsigned ldA, double const* B, unsigned ldB, double beta, double* C, unsigned ldC) {
  for (unsigned j = 0; j < N; ++j) {
    for (unsigned i = 0; i < M; ++i) {
      double cij = 0.0;
      for (unsigned k = 0; k < K; ++k) {
        cij += A[k*ldA + i] * B[j*ldB + k];
      }
      C[j*ldC + i] = alpha * cij + beta * C[j*ldC + i];
    }
  }
}
// same as DGEMM() for Fallback
void opt_DGEMM2(unsigned M, unsigned N, unsigned K, double alpha, double const* A, unsigned ldA, double const* B, unsigned ldB, double beta, double* C, unsigned ldC) {
  for (unsigned j = 0; j < N; ++j) {
    for (unsigned i = 0; i < M; ++i) {
      double cij = 0.0;
      for (unsigned k = 0; k < K; ++k) {
        cij += A[k*ldA + i] * B[j*ldB + k];
      }
      C[j*ldC + i] = alpha * cij + beta * C[j*ldC + i];
    }
  }
}
// same as DGEMM() for Fallback
void opt_DGEMM4A(unsigned M, unsigned N, unsigned K, double alpha, double const* A, unsigned ldA, double const* B, unsigned ldB, double beta, double* C, unsigned ldC) {
  for (unsigned j = 0; j < N; ++j) {
    for (unsigned i = 0; i < M; ++i) {
      double cij = 0.0;
      for (unsigned k = 0; k < K; ++k) {
        cij += A[k*ldA + i] * B[j*ldB + k];
      }
      C[j*ldC + i] = alpha * cij + beta * C[j*ldC + i];
    }
  }
}
// same as DGEMM() for Fallback
void opt_DGEMM4B(unsigned M, unsigned N, unsigned K, double alpha, double const* A, unsigned ldA, double const* B, unsigned ldB, double beta, double* C, unsigned ldC) {
  for (unsigned j = 0; j < N; ++j) {
    for (unsigned i = 0; i < M; ++i) {
      double cij = 0.0;
      for (unsigned k = 0; k < K; ++k) {
        cij += A[k*ldA + i] * B[j*ldB + k];
      }
      C[j*ldC + i] = alpha * cij + beta * C[j*ldC + i];
    }
  }
}

#else

// optimized

void opt_DGEMM1(unsigned M, unsigned N, unsigned K, double alpha, double const* A, unsigned ldA, double const* B, unsigned ldB, double beta, double* C, unsigned ldC) {

  __m512d A_reg;
  __m512d B_reg[PADDED_GLOBAL_MATRICES_LENGTH/8];
  __m512d C_reg;

  for (unsigned j = 0; j < N; ++j) {

    for (unsigned k = 0; k < PADDED_GLOBAL_MATRICES_LENGTH/8; ++k) {
      B_reg[k] = _mm512_load_pd(&B[j*ldB + k*8]);
    }

    for (unsigned i = 0; i < M; ++i) {
      double cij = 0.0;
      C_reg = _mm512_set1_pd(cij);
      for (unsigned k = 0; k < PADDED_GLOBAL_MATRICES_LENGTH/8; ++k) {
        A_reg = _mm512_load_pd(&A[i*PADDED_GLOBAL_MATRICES_LENGTH + k*8]);
        C_reg = _mm512_fmadd_pd(A_reg, B_reg[k], C_reg);
      }
      cij = _mm512_reduce_add_pd(C_reg);
      
      // beta here always 0
      C[i*8 + j] = alpha * cij;
    }
  }
}

void opt_DGEMM2(unsigned M, unsigned N, unsigned K, double alpha, double const* A, unsigned ldA, double const* B, unsigned ldB, double beta, double* C, unsigned ldC) {

  __m512d A_reg;
  __m512d B_reg;
  __m512d C_reg;
  double cij;

  for (unsigned j = 0; j < N; ++j) {
    B_reg = _mm512_load_pd(&B[j*ldB]);
    for (unsigned i = 0; i < M; ++i) {
      
      A_reg = _mm512_load_pd(&A[i*8]);
      C_reg = _mm512_mul_pd(A_reg, B_reg);
      cij = _mm512_reduce_add_pd(C_reg);

      // beta here always 1
      C[j*ldC + i] = alpha * cij + C[j*ldC + i];
    }
  }
}

void opt_DGEMM4A(unsigned M, unsigned N, unsigned K, double alpha, double const* A, unsigned ldA, double const* B, unsigned ldB, double beta, double* C, unsigned ldC) {

  // B only with values in 0/1 and 1/0

  __m512d A_reg;
  __m512d B_reg;
  __m512d C_reg;
  double cij;

  // j = 0
    B_reg = _mm512_load_pd(&B[0]);
    for (unsigned i = 0; i < M; ++i) {
      
      A_reg = _mm512_load_pd(&A[i*8]);
      C_reg = _mm512_mul_pd(A_reg, B_reg);
      cij = _mm512_reduce_add_pd(C_reg);

      // beta here always 1
      C[i] = alpha * cij + C[i];
    }
  // j = 1
    B_reg = _mm512_load_pd(&B[ldB]);
    for (unsigned i = 0; i < M; ++i) {
      
      A_reg = _mm512_load_pd(&A[i*8]);
      C_reg = _mm512_mul_pd(A_reg, B_reg);
      cij = _mm512_reduce_add_pd(C_reg);

      // beta here always 1
      C[ldC + i] = alpha * cij + C[ldC + i];
    }

}

void opt_DGEMM4B(unsigned M, unsigned N, unsigned K, double alpha, double const* A, unsigned ldA, double const* B, unsigned ldB, double beta, double* C, unsigned ldC) {

  // B has only values in 0/2 and 2/0

  __m512d A_reg;
  __m512d B_reg;
  __m512d C_reg;
  double cij;

  // j = 0
    B_reg = _mm512_load_pd(&B[0]);
    for (unsigned i = 0; i < M; ++i) {
      
      A_reg = _mm512_load_pd(&A[i*8]);
      C_reg = _mm512_mul_pd(A_reg, B_reg);
      cij = _mm512_reduce_add_pd(C_reg);

      // beta here always 1
      C[i] = alpha * cij + C[i];
    }
  // j = 2
    B_reg = _mm512_load_pd(&B[2*ldB]);
    for (unsigned i = 0; i < M; ++i) {
      
      A_reg = _mm512_load_pd(&A[i*8]);
      C_reg = _mm512_mul_pd(A_reg, B_reg);
      cij = _mm512_reduce_add_pd(C_reg);

      // beta here always 1
      C[2*ldC + i] = alpha * cij + C[2*ldC + i];
    }

}

#endif

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <limits>

#include "Stopwatch.h"

void dgemm(double* A, double* B, double* C) {
  for (int n = 0; n < N; ++n) {
    for (int m = 0; m < M; ++m) {
      for (int k = 0; k < K; ++k) {
        C[n*M + m] += A[k*M + m] * B[n*K + k];
      }
    }
  }
}

void dgemm_opt(double* A, double* B, double* C) {
  
  /** =================
   *         TODO
   *  ================= */
  
}

int main(int argc, char** argv) {
  int repetitions = 10000;
  if (argc > 1) {
    repetitions = atoi(argv[1]);
  }
  
  /** Allocate memory */
  double* A, *B, *C, *A_test, *B_test, *C_test;
  
  posix_memalign(reinterpret_cast<void**>(&A),      ALIGNMENT, M*K*sizeof(double));
  posix_memalign(reinterpret_cast<void**>(&B),      ALIGNMENT, K*N*sizeof(double));
  posix_memalign(reinterpret_cast<void**>(&C),      ALIGNMENT, M*N*sizeof(double));
  posix_memalign(reinterpret_cast<void**>(&A_test), ALIGNMENT, M*K*sizeof(double));
  posix_memalign(reinterpret_cast<void**>(&B_test), ALIGNMENT, K*N*sizeof(double));
  posix_memalign(reinterpret_cast<void**>(&C_test), ALIGNMENT, M*N*sizeof(double));

  for (int j = 0; j < K; ++j) {
    for (int i = 0; i < M; ++i) {
      A[j*M + i] = i + j;
    }
  }
  for (int j = 0; j < N; ++j) {
    for (int i = 0; i < K; ++i) {
      B[j*K + i] = (K-i) + (N-j);
    }
  }
  memset(C, 0, M*N*sizeof(double));
  memcpy(A_test, A, M*K*sizeof(double));
  memcpy(B_test, B, K*N*sizeof(double));
  memset(C_test, 0, M*N*sizeof(double));
  
  /** Check correctness of optimised dgemm */
  #pragma noinline
  {
    dgemm(A, B, C);
    dgemm_opt(A_test, B_test, C_test);
  }

  double error = 0.0;
  for (int i = 0; i < M*N; ++i) {
    double diff = C[i] - C_test[i];
    error += diff*diff;
  }
  error = sqrt(error);
  if (error > std::numeric_limits<double>::epsilon()) {
    printf("Optimised DGEMM is incorrect. Error: %e\n", error);
    return -1;
  }
  
  /** Test performance of optimised dgemm */
  
  #pragma noinline
  dgemm_opt(A, B, C);
  
  Stopwatch stopwatch;
  stopwatch.start();

  #pragma noinline
  for (int r = 0; r < repetitions; ++r) {
    dgemm_opt(A, B, C);
    __asm__ __volatile__("");
  }
  
  double time = stopwatch.stop();
  printf("%lf ms, %lf GFLOP/s\n", time * 1.0e3, repetitions*2.0*M*N*K/time * 1.0e-9);
  
  /** Clean up */
  
  free(A); free(B); free(C);
  free(A_test); free(B_test); free(C_test);

  return 0;
}

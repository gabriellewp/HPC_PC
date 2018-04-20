#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <limits>
#include <algorithm>

#include <immintrin.h>
#include <omp.h>

#include "Stopwatch.h"

/** Syntax: M x N, ld S
 *  A M x N sub-block of of S x something matrix
 *  in column-major layout.
 *  That is C_ij := C[j*S + i], i=0,...,M-1,  j=0,...,N-1
 */


void dgemm(double* A, double* B, double* C, int S) {
  for (int n = 0; n < S; ++n) {
    for (int m = 0; m < S; ++m) {
      for (int k = 0; k < S; ++k) {
        C[n*S + m] += A[k*S + m] * B[n*S + k];
      }
    }
  }
}

/**
 * A: M x K, ld M
 * B: K x N, ld K
 * C: M x N, ld S
 */
void microkernel(double* A, double* B, double* C, int S) {
	
	/** =================
	*         TODO
	*  ================= */ 
	
	//printf("CALL microkernel\n");
	
	//              !! M durch 8 teilbar !!
	// beachten, sonst wird A nicht komplett durchgegangen
	
	__m512d A_reg;
	__m512d B_reg;
	
	__m512d C_reg[M*N/8];
	
	for(int i = 0; i < N; ++i){
		for(int j = 0; j < M/8; ++j){
			C_reg[(M/8)*i+j] = _mm512_load_pd(&C[i*S+j*8]);
		}
	}
	
	for(int i = 0; i < K; ++i){
		for(int j = 0; j < M/8; ++j){
			A_reg = _mm512_load_pd(&A[i*MC+j*8]);
			for(int n = 0; n < N; ++n){
				B_reg = _mm512_set1_pd(B[n*K+i]);
				C_reg[n*(M/8)+j] = _mm512_fmadd_pd(A_reg, B_reg, C_reg[n*(M/8)+j]);
			}
		}
	}
	
	for(int i = 0; i < N; ++i){
		for(int j = 0; j < M/8; ++j){
			_mm512_store_pd(&C[i*S+j*8], C_reg[(M/8)*i+j]);
		}
	}
	
}

/**
 * A: MC x K, ld S
 * B:  K x S, ld K
 * C: MC x S, ld S
 */
void GEBP(double* A, double* B, double* C, double* A_pack, int S, int threadsPerTeam){
	
	/** =================
	*         TODO
	*  ================= */
	
	//printf("CALL GEBP\n");
	
	//                 !! MC/M % threadsPerTeam == 0 !!
	// beachten, sonst wird nicht alles komplett oder zuviel durchgegangen
	
	#pragma omp parallel num_threads(threadsPerTeam)
	{
		int rank = omp_get_thread_num();
		int MsPerThread = (MC/M)/threadsPerTeam;
		
		for(int i = 0; i < MsPerThread; ++i){
			
			int part = MsPerThread*rank+i;
			
			for(int k = 0; k < S/N; ++k){
				microkernel(&A_pack[part*M], &B[k*K*N], &C[k*S*N+part*M], S);
			}
		}
	}
	
}

/**
 * A: S x K, ld S
 * B: K x S, ld S
 * C: S x S, ld S
 */
void GEPP(double* A, double* B, double* C, double** A_pack, double* B_pack, int S, int nTeams, int threadsPerTeam){
	
	/** =================
	*         TODO
	*  ================= */
	
	//printf("CALL GEPP\n");
	
	//               !! nTeams == S/MC !!
	// beachten, sonst wird A nicht komplett durchgegangen
	
	//for (int t = 0; t < nTeams; ++t)
	#pragma omp parallel num_threads(nTeams)
	{
		int t = omp_get_thread_num();
		
		for(int i = 0; i < K; ++i){
			for(int j = 0; j < MC; ++j){
				A_pack[t][i*MC + j] = A[S*i+t*MC+j];
			}
		}
		
		GEBP(A, B_pack, &C[t*MC], A_pack[t], S, threadsPerTeam);
	}
	
}

/**
 * A: S x S, ld S
 * B: S x S, ld S
 * C: S x S, ld S
 */
void GEMM(double* A, double* B, double* C, double** A_pack, double* B_pack, int S, int nTeams = 1, int threadsPerTeam = 1) {
	
	/** =================
	*         TODO
	*  ================= */
	
	//printf("CALL GEMM\n");
	
	for(int i = 0; i < S/K; ++i){
		
		for(int j = 0; j < S; ++j){
			for(int k = 0; k < K; ++k){
				B_pack[j*K + k] = B[i*K+j*S+k];
			}
		}
		
		GEPP(&A[i*S*K], B_pack, C, A_pack, B_pack, S, nTeams, threadsPerTeam);
	}

}

int main(int argc, char** argv) {
  int S = 4096;
  bool test = true;
  int threadsPerTeam = 1;
  int nRepeat = 10;
  if (argc <= 1) {
    printf("Usage: dgemm <S> <test> <threadsPerTeam> <repetitions>");
    return -1;
  }
  if (argc > 1) {
    S = atoi(argv[1]);
  }
  if (argc > 2) {
    test = atoi(argv[2]) != 0;
  }
  if (argc > 3) {
    threadsPerTeam = atoi(argv[3]);
  }
  if (argc > 4) {
    nRepeat = atoi(argv[4]);
  }

  omp_set_nested(1);

  int nThreads, nTeams;
  #pragma omp parallel
  #pragma omp master
  {
    nThreads = omp_get_num_threads(); 
  }
  threadsPerTeam = std::min(threadsPerTeam, nThreads);
  nTeams = nThreads / threadsPerTeam;
  
  printf("Teams: %d, Thread-per-Team: %d\n", nTeams, threadsPerTeam);
  
  /** Allocate memory */
  double* A, *B, *C, *A_test, *B_test, *C_test, **A_pack, *B_pack, *C_aux;
  
  posix_memalign(reinterpret_cast<void**>(&A),      ALIGNMENT, S*S*sizeof(double));
  posix_memalign(reinterpret_cast<void**>(&B),      ALIGNMENT, S*S*sizeof(double));
  posix_memalign(reinterpret_cast<void**>(&C),      ALIGNMENT, S*S*sizeof(double));
  posix_memalign(reinterpret_cast<void**>(&A_test), ALIGNMENT, S*S*sizeof(double));
  posix_memalign(reinterpret_cast<void**>(&B_test), ALIGNMENT, S*S*sizeof(double));
  posix_memalign(reinterpret_cast<void**>(&C_test), ALIGNMENT, S*S*sizeof(double));

  posix_memalign(reinterpret_cast<void**>(&A_pack), ALIGNMENT, nTeams*sizeof(double*));
  for (int t = 0; t < nTeams; ++t) {
    posix_memalign(reinterpret_cast<void**>(&A_pack[t]), ALIGNMENT, MC*K*sizeof(double));
  }
  posix_memalign(reinterpret_cast<void**>(&B_pack), ALIGNMENT,  K*S*sizeof(double));

  #pragma omp parallel for
  for (int j = 0; j < S; ++j) {
    for (int i = 0; i < S; ++i) {
      A[j*S + i] = i + j;
      B[j*S + i] = (S-i) + (S-j);
      C[j*S + i] = 0.0;
    }
  }
  memcpy(A_test, A, S*S*sizeof(double));
  memcpy(B_test, B, S*S*sizeof(double));
  memset(C_test, 0, S*S*sizeof(double));
  
  /** Check correctness of optimised dgemm */
  if (test) {
    #pragma noinline
    {
      dgemm(A_test, B_test, C_test, S);
      GEMM(A, B, C, A_pack, B_pack, S, nTeams, threadsPerTeam);
    }

    double error = 0.0;
    for (int i = 0; i < S*S; ++i) {
      double diff = C[i] - C_test[i];
	  //printf("%f -- %f", C[i], C_test[i]);
      error += diff*diff;
    }
    error = sqrt(error);
    if (error > std::numeric_limits<double>::epsilon()) {
      printf("Optimised DGEMM is incorrect. Error: %e\n", error);
      return -1;
    }else{
      printf("Optimised DGEMM is CORRECT :)\n");
	}
  }

  Stopwatch stopwatch;
  double time;

  /** Test performance of microkernel */

  stopwatch.start();
  for (int i = 0; i < 10000; ++i) {
    #pragma noinline
    microkernel(A, B, C, S);
    __asm__ __volatile__("");
  }
  time = stopwatch.stop();
  printf("Microkernel: %lf ms, %lf GFLOP/s\n", time*1.0e3, 10000*2.0*M*N*K/time * 1.0e-9);

  /** Test performance of GEBP */

  stopwatch.start();
  for (int i = 0; i < nRepeat; ++i) {
    #pragma noinline
    GEBP(A, B, C, A_pack[0], S, threadsPerTeam);
    __asm__ __volatile__("");
  }
  time = stopwatch.stop();
  printf("GEBP: %lf ms, %lf GFLOP/s\n", time*1.0e3, nRepeat*2.0*MC*S*K/time * 1.0e-9);

  /** Test performance of optimised GEMM */

  stopwatch.start();
  for (int i = 0; i < nRepeat; ++i) {
    #pragma noinline
    GEMM(A, B, C, A_pack, B_pack, S, nTeams, threadsPerTeam);
    __asm__ __volatile__("");
  }  
  time = stopwatch.stop();
  printf("GEMM: %lf ms, %lf GFLOP/s\n", time * 1.0e3, nRepeat*2.0*S*S*S/time * 1.0e-9);
  
  /** Clean up */
  
  free(A); free(B); free(C);
  free(A_test); free(B_test); free(C_test);
  for (int t = 0; t < nTeams; ++t) {
    free(A_pack[t]);
  }
  free(A_pack); free(B_pack);

  return 0;
}

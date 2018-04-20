#include <cstdlib>
#include <omp.h>
#ifdef LIKWID_PERFMON
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif
#define N 10000
#define R 10

void dtrmv(double const* A, double const* x, double* y) {
  #pragma omp parallel
  {
     LIKWID_MARKER_START("compute");
        #pragma omp for
        for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                 y[i] += A[i*N + j] * x[j];
                }
         }

    LIKWID_MARKER_STOP("compute");
  }
}

int main() {
  double *A, *x, *y;
  LIKWID_MARKER_INIT;
  #pragma omp parallel
  {
    LIKWID_MARKER_THREADINIT;
  }
  posix_memalign(reinterpret_cast<void**>(&A), 64, N*N*sizeof(double));
  posix_memalign(reinterpret_cast<void**>(&x), 64,   N*sizeof(double));
  posix_memalign(reinterpret_cast<void**>(&y), 64,   N*sizeof(double));

  for (unsigned i = 0; i < N; ++i) {
    for (unsigned j = 0; j < N; ++j) {
      A[i*N + j] = i + j;
    }
    x[i] = i;
  }

  for (unsigned r = 0; r < R; ++r) {
    dtrmv(A, x, y);
  }

  free(A); free(x); free(y);
  LIKWID_MARKER_CLOSE;
  return 0;
}

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "Stopwatch.h"

double function(double x) {
  return x + sin(x);
}

void kernel1(double* in, double* out, int N) {
  for (unsigned i = 0; i < N; ++i) {
    out[i] = function(in[i]);
  }
}

void kernel2(double* in, double* out, int N) {
  for (unsigned i = 0; i < N; ++i) {
    out[i] = 3.14 * in[i] - 5.0;
  }
}


int main() {
  int N = 1000000;
  int R = 10000;
  
  double *in, *out;
  posix_memalign(reinterpret_cast<void**>( &in), 64, N*sizeof(double));
  posix_memalign(reinterpret_cast<void**>(&out), 64, N*sizeof(double));
  
  Stopwatch sw;
  sw.start();

  for (unsigned r = 0; r < R; ++r) {
    if (drand48() < 0.2) {
      kernel1(in, out, N);
    } else {
      kernel2(in, out, N);
    }
  }

  double time = sw.stop();
  printf("Elapsed time: %lf\n", time);
  
  free(in);
  free(out);
  
  return 0;
}


#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include "Stopwatch.h"
int main(int argc, char **argv) {

  long int iterations = 100000;
  double stepsize = 0;
  double x = 0;
  double y_step = 0;
  double pi_4 = 0;
  double end_result = 0;
  stepsize = (double)1/iterations;
  long int step = 0;
  Stopwatch stopwatch;
  stopwatch.start();
  
 for(step = 0; step < iterations; step++) {
   x = (double) step/iterations;
   y_step = (double) 1/(1+(x*x));
   pi_4 =  pi_4 + (double) (y_step*stepsize);
}

  end_result = 4*stepsize*pi_4;
  printf("%.6f\n", end_result);

  double time = stopwatch.stop();
  printf("%lf ms\n", time );

  return 0;
}

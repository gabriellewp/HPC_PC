#include <cstdlib>
#include <cstdio>
#include "Stopwatch.h"
#include <likwid.h>


class RigidBody {
public:
  double x,y,z;

  RigidBody() {
    x = drand48();
    y = drand48();
    z = drand48();
  }

  inline void move(double dx, double dy, double dz) {
    x += dx;
    y += dy;
    z += dz;
  }
};

void gravity(double dt, RigidBody** bodies, int N) {
  LIKWID_MARKER_START("Compute");
  for (int n = 0; n < N; ++n) {
    bodies[n]->move(0.0, 0.0, 0.5 * 9.81 * dt * dt);
  }
  LIKWID_MARKER_STOP("Compute");
}

int main() {
  int N = 10000;
  int T = 10;
  LIKWID_MARKER_INIT;

  RigidBody** bodies = new RigidBody*[N];

  for (int i = 0; i < N; ++i) {
    bodies[i] = new RigidBody[i];
  }

  Stopwatch sw;
  sw.start();

  for (int t = 0; t < T; ++t) {
    gravity(0.001, bodies, N);
  }

  double time = sw.stop();
  printf("Time: %lf s, BW: %lf\n MB/s", time, T*N*sizeof(double)*1e-6 / time);

  for (int i = 0; i < N; ++i) {
    delete bodies[i];
  }

  delete[] bodies;

  LIKWID_MARKER_CLOSE;
  return 0;
}

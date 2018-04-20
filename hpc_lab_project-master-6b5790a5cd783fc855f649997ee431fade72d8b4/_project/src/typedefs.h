#ifndef TYPEDEFS_H_
#define TYPEDEFS_H_

#include <cmath>
#include "constants.h"

typedef double DegreesOfFreedom[NUMBER_OF_DOFS];

struct GlobalConstants {
  double hx;
  double hy;
  int X;
  int Y;
  double maxTimestep;
  double endTime;
  
  double c1; // -globals.hx / (globals.hx * globals.hy)
  double c2; // -globals.hy / (globals.hx * globals.hy)
  double c3; // 1. / (globals.hx*globals.hy)
  double c4; // -1.0 / globals.hx
  double c5; // -1.0 / globals.hy
  double c6; // 1.0 / globals.hx
  double c7; // 1.0 / globals.hy
};

struct Material {
  double K0;
  double rho0;
  
  inline double wavespeed() const { return sqrt(K0/rho0); }
};

struct SourceTerm {
  SourceTerm() : x(-1), y(-1) {} // -1 == invalid
  int x;
  int y;
  double phi[NUMBER_OF_BASIS_FUNCTIONS];
  double (*antiderivative)(double);
  int quantity;
};

#endif // TYPEDEFS_H_

#include "Simulator.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include "omp.h"
#include "constants.h"
#include "Kernels.h"
#include "Model.h"
#include "GlobalMatrices.h"
#include "mpi.h"

double determineTimestep(double hx, double hy, Grid<Material>& materialGrid)
{
  double maxWaveSpeed = 0.0;
  for (int y = 0; y < materialGrid.Y(); ++y) {
    for (int x = 0; x < materialGrid.X(); ++x) {
      maxWaveSpeed = std::max(maxWaveSpeed, materialGrid.get(x, y).wavespeed());
    }
  }
  
  return 0.25 * std::min(hx, hy)/((2*CONVERGENCE_ORDER-1) * maxWaveSpeed);
}

int simulate( GlobalConstants const&  globals,
              Grid<Material>&         materialGrid,
              Grid<DegreesOfFreedom>& degreesOfFreedomGrid,
              WaveFieldWriter&        waveFieldWriter,
              SourceTerm&             sourceterm  )
{
  
  Grid<DegreesOfFreedom> timeIntegratedGrid(globals.X, globals.Y);
  
  double time;
  int step = 0;
  for (time = 0.0; time < globals.endTime; time += globals.maxTimestep) {
    waveFieldWriter.writeTimestep(time, degreesOfFreedomGrid);
  
    double timestep = std::min(globals.maxTimestep, globals.endTime - time);
    
    #pragma omp parallel for collapse(2)
    for (int y = 0; y < globals.Y; ++y) {
      for (int x = 0; x < globals.X; ++x) {
        double Aplus[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES];
        double rotatedAplus[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES];
        
        Material& material = materialGrid.get(x, y);
        DegreesOfFreedom& degreesOfFreedom = degreesOfFreedomGrid.get(x, y);
        DegreesOfFreedom& timeIntegrated = timeIntegratedGrid.get(x, y);
        
        computeAder(timestep, globals, material, degreesOfFreedom, timeIntegrated);
        
        computeVolumeIntegral(globals, material, timeIntegrated, degreesOfFreedom);

        computeAplus(material, materialGrid.get(x, y-1), Aplus);
        rotateFluxSolver(0., -1., Aplus, rotatedAplus);
        computeFlux(globals.c1, GlobalMatrices::Fxm0, rotatedAplus, timeIntegrated, degreesOfFreedom);
        
        computeAplus(material, materialGrid.get(x, y+1), Aplus);
        rotateFluxSolver(0., 1., Aplus, rotatedAplus);
        computeFlux(globals.c1, GlobalMatrices::Fxm1, rotatedAplus, timeIntegrated, degreesOfFreedom);
        
        computeAplus(material, materialGrid.get(x-1, y), Aplus);
        rotateFluxSolver(-1., 0., Aplus, rotatedAplus);
        computeFlux(globals.c2, GlobalMatrices::Fym0, rotatedAplus, timeIntegrated, degreesOfFreedom);
        
        computeAplus(material, materialGrid.get(x+1, y), Aplus);
        rotateFluxSolver(1., 0., Aplus, rotatedAplus);
        computeFlux(globals.c2, GlobalMatrices::Fym1, rotatedAplus, timeIntegrated, degreesOfFreedom);
      }
    }
    #pragma omp parallel for collapse(2)
    for (int y = 0; y < globals.Y; ++y) {
      for (int x = 0; x < globals.X; ++x) {
        double Aplus[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES];
        double rotatedAplus[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES];

        Material& material = materialGrid.get(x, y);
        DegreesOfFreedom& degreesOfFreedom = degreesOfFreedomGrid.get(x, y);

        computeAminus(material, materialGrid.get(x, y-1), Aplus);
        rotateFluxSolver(0., -1., Aplus, rotatedAplus);
        computeFlux(globals.c1, GlobalMatrices::Fxp0, rotatedAplus, timeIntegratedGrid.get(x, y-1), degreesOfFreedom);
        
        computeAminus(material, materialGrid.get(x, y+1), Aplus);
        rotateFluxSolver(0., 1., Aplus, rotatedAplus);
        computeFlux(globals.c1, GlobalMatrices::Fxp1, rotatedAplus, timeIntegratedGrid.get(x, y+1), degreesOfFreedom);
        
        computeAminus(material, materialGrid.get(x-1, y), Aplus);
        rotateFluxSolver(-1., 0., Aplus, rotatedAplus);
        computeFlux(globals.c2, GlobalMatrices::Fyp0, rotatedAplus, timeIntegratedGrid.get(x-1, y), degreesOfFreedom);
        
        computeAminus(material, materialGrid.get(x+1, y), Aplus);
        rotateFluxSolver(1., 0., Aplus, rotatedAplus);
        computeFlux(globals.c2, GlobalMatrices::Fyp1, rotatedAplus, timeIntegratedGrid.get(x+1, y), degreesOfFreedom);
      }
    }
    
    if (sourceterm.x >= 0 && sourceterm.x < globals.X && sourceterm.y >= 0 && sourceterm.y < globals.Y) {
      DegreesOfFreedom& degreesOfFreedom = degreesOfFreedomGrid.get(sourceterm.x, sourceterm.y);
      double timeIntegral = (*sourceterm.antiderivative)(time + timestep) - (*sourceterm.antiderivative)(time);
      for (unsigned b = 0; b < NUMBER_OF_BASIS_FUNCTIONS; ++b) {
        degreesOfFreedom[sourceterm.quantity * NUMBER_OF_BASIS_FUNCTIONS + b] += globals.c3 * timeIntegral * sourceterm.phi[b];
      }
    }
    
    ++step;
    if (step % 100 == 0) {
      std::cout << "At time / timestep: " << time << " / " << step << std::endl;
    }
  }
  
  waveFieldWriter.writeTimestep(globals.endTime, degreesOfFreedomGrid, true);
  return step;
}

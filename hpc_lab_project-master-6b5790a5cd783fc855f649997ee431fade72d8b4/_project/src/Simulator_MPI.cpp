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
  int ntasks, rank, left, right;
  
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Request request[4];
  MPI_Status status[4];

  int flag[4] = {0, 0, 0, 0};
  
  /*new solution into old*/
  left = rank - 1;
  right = rank + 1;
    if(rank == 0) left = ntasks-1;
    if(rank == ntasks-1) right = 0;
  
   
  int stripsize=( globals.X)/ ntasks;
  int mlo, mhi;
  Grid<DegreesOfFreedom> timeIntegratedGrid(globals.X, globals.Y);
  Grid<DegreesOfFreedom> timeIntegratedGridBuffLeft(1, globals.Y);
  Grid<DegreesOfFreedom> timeIntegratedGridBuffRight(1, globals.Y);
    //define the halo index
    if(rank ==0){
        mlo =0;
        mhi = mlo + stripsize-1;
    }else if(rank==ntasks-1){
        mlo =(rank * stripsize) -1;
        mhi = globals.X -1;
    }else{
        mlo = (rank * stripsize)-1;
        mhi = mlo + stripsize-1;
    }

  double time;
  int step = 0;
  for (time = 0.0; time < globals.endTime; time += globals.maxTimestep) {
    for(int iflag=0;iflag<4;iflag++){
        flag[iflag]=0;
    }
    waveFieldWriter.writeTimestep(time, degreesOfFreedomGrid); 
    double timestep = std::min(globals.maxTimestep, globals.endTime - time);
   
    //local computation, no communication
    #pragma omp parallel for collapse(2)
    for (int y = 0; y < globals.Y; ++y) {
      for (int x = mlo; x <= mhi; ++x) {
        double Aplus[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES];
        double rotatedAplus[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES];
        
        Material& material = materialGrid.get(x, y);
        DegreesOfFreedom& degreesOfFreedom = degreesOfFreedomGrid.get(x, y);
        DegreesOfFreedom& timeIntegrated = timeIntegratedGrid.get(x, y);
        

        memcpy(&timeIntegratedGridBuffLeft.get(0,0), &timeIntegratedGrid.get(mlo,0), globals.Y);
        memcpy(&timeIntegratedGridBuffRight.get(0,0), &timeIntegratedGrid.get(mhi,0), globals.Y);      
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

    //entering second loop, start to communicate with neighbor
    
     //receive from the left
    MPI_Irecv(&timeIntegratedGridBuffLeft.get(0,0), globals.Y, MPI_DOUBLE, left, 1, MPI_COMM_WORLD, &request[0]);
    //receive from the right
    MPI_Irecv(&timeIntegratedGridBuffRight.get(0,0), globals.Y, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &request[1]);
    
    //compute left and right halo column
    int i = 2, x=0;
    while(i < 2){
       if(i==0){
           x = mlo;
        }else{
           x = mhi;
        }
 
    for (int y = 0; y < globals.Y; ++y) {
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
        computeFlux(globals.c2, GlobalMatrices::Fyp0, rotatedAplus, timeIntegratedGridBuffLeft.get(0, y), degreesOfFreedom);
        memcpy(&timeIntegratedGrid.get(x-1,y),&timeIntegratedGridBuffLeft.get(0,y),1);     
    
        computeAminus(material, materialGrid.get(x+1, y), Aplus);
        rotateFluxSolver(1., 0., Aplus, rotatedAplus);
        computeFlux(globals.c2, GlobalMatrices::Fyp1, rotatedAplus, timeIntegratedGridBuffRight.get(0, y), degreesOfFreedom);
        memcpy(&timeIntegratedGrid.get(x+1,y),&timeIntegratedGridBuffRight.get(0,y),1);
    }
     
     i++;
    }
    
    //make sure all receive finished
      MPI_Wait(&request[0], &status[0]);
      MPI_Wait(&request[1], &status[1]);
      MPI_Test(&request[0], &flag[0], &status[0]);
      MPI_Test(&request[1], &flag[1], &status[1]);

    //copy the boundary result into the buffLeft and buffright
    memcpy(&timeIntegratedGridBuffLeft.get(0,0), &timeIntegratedGrid.get(mlo,0), globals.Y);
    memcpy(&timeIntegratedGridBuffRight.get(0,0), &timeIntegratedGrid.get(mhi,0), globals.Y);
    

    //send to neighbour 
    //send to the right
    MPI_Isend(&timeIntegratedGridBuffRight.get(0,0), globals.Y, MPI_DOUBLE, right, 1, MPI_COMM_WORLD, &request[2]);
    //send to the left
    MPI_Isend(&timeIntegratedGridBuffLeft.get(0,0), globals.Y, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &request[3]);
    
    //compute the rest of inner domain
    #pragma omp parallel for collapse(2)
    for (int y = 0; y < globals.Y; ++y) {
      for (int x = mlo+1; x < mhi; ++x) { 
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
    MPI_Wait(&request[2], &status[2]);
    MPI_Wait(&request[3], &status[3]);

    MPI_Test(&request[2], &flag[2], &status[2]);
    MPI_Test(&request[3], &flag[3], &status[3]);   



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
  for(int i =0; i <4 ; i++){
       printf("flag %d", flag[i]);

  }

  waveFieldWriter.writeTimestep(globals.endTime, degreesOfFreedomGrid, true);

  return step;
}

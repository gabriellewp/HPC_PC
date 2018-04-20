#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void my_bcast(void* data, int count, MPI_Datatype datatype, int root,
              MPI_Comm communicator) {
  int world_rank;
  MPI_Comm_rank(communicator, &world_rank);
  int world_size;
  MPI_Comm_size(communicator, &world_size);

  if (world_rank == root) {
    // If we are the root process, send our data to everyone
    int i;
    for (i = 0; i < world_size; i++) {
      if (i != world_rank) {
        MPI_Send(data, count, datatype, i, 0, communicator);
      }
    }
  } else {
    MPI_Recv(data, count, datatype, root, 0, communicator, MPI_STATUS_IGNORE);
  }
}

int main(int argc, char** argv) {
  MPI_Init(NULL, NULL);
  int count = atoi(argv[1]);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);


  double my_bcast_time = 0.0;
  double mpibcast_time = 0.0;
  double buf[count];
  int i;
   for (i = 0; i < count; i++)
           buf[i] = 0.0;

  MPI_Barrier(MPI_COMM_WORLD);
  my_bcast_time -= MPI_Wtime();
  my_bcast(&buf, count, MPI_INT, 0, MPI_COMM_WORLD);
  // Synchronize again before obtaining final time
  MPI_Barrier(MPI_COMM_WORLD);
  my_bcast_time += MPI_Wtime();

  // Time MPI_Bcast
  MPI_Barrier(MPI_COMM_WORLD);
  mpibcast_time -= MPI_Wtime();
  MPI_Bcast(&buf, count, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  mpibcast_time += MPI_Wtime();

  if (world_rank == 0) {
    printf("mybcast time = %lf\n", my_bcast_time);
    printf("MPI_Bcast time = %lf\n", mpibcast_time);
  }



  MPI_Finalize();
}

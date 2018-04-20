/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the training material of the master's course           *
* Scientific Computing                                                        *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <x86intrin.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <sys/time.h>

#include <mpi.h>

#include <omp.h>

MPI_Datatype myVectorTypeGridWide;
MPI_Datatype myVectorTypeGrid;
MPI_Datatype myVectorType1Col;

/// store number of grid points in one dimension
std::size_t grid_points_1d = 0;
std::size_t grid_points_x = 0;
std::size_t grid_points_y = 0;

int rank = 0;
int size = 0;

/// store begin timestep
struct timeval begin;
/// store end timestep
struct timeval end;

/**
 * initialize and start timer
 */
void timer_start()
{
	gettimeofday(&begin,(struct timezone *)0);
}

/**
 * stop timer and return measured time
 *
 * @return measured time
 */
double timer_stop()
{
	gettimeofday(&end,(struct timezone *)0);
	double seconds, useconds;
	double ret, tmp;

	if (end.tv_usec >= begin.tv_usec)
	{
		seconds = (double)end.tv_sec - (double)begin.tv_sec;
		useconds = (double)end.tv_usec - (double)begin.tv_usec;
	}
	else
	{
		seconds = (double)end.tv_sec - (double)begin.tv_sec;
		seconds -= 1;					// Correction
		useconds = (double)end.tv_usec - (double)begin.tv_usec;
		useconds += 1000000;			// Correction
	}

	// get time in seconds
	tmp = (double)useconds;
	ret = (double)seconds;
	tmp /= 1000000;
	ret += tmp;

	return ret;
}

/**
 * stores a given grid into a file
 * 
 * @param grid the grid that should be stored
 * @param filename the filename
 */
void store_grid(double* grid, std::string filename)
{
	std::fstream filestr;
	filestr.open (filename.c_str(), std::fstream::out);
	
	// calculate mesh width 
	double mesh_width = 1.0/((double)(grid_points_1d-1));

	// store grid incl. boundary points
	for (int i = 0; i < grid_points_1d; i++)
	{
		for (int j = 0; j < grid_points_1d; j++)
		{
			filestr << mesh_width*i << " " << mesh_width*j << " " << grid[(i*grid_points_1d)+j] << std::endl;
		}
		
		filestr << std::endl;
	}

	filestr.close();
}

/**
 * calculate the grid's initial values for given grid points
 *
 * @param x the x-coordinate of a given grid point
 * @param y the y-coordinate of a given grid point
 *
 * @return the initial value at position (x,y)
 */
double eval_init_func(double x, double y)
{
	return (x*x)*(y*y);
}

/**
 * initializes a given grid: inner points are set to zero
 * boundary points are initialized by calling eval_init_func
 *
 * @param grid the grid to be initialized
 */
void init_grid(double* grid)
{
	// set all points to zero
	for (int i = 0; i < grid_points_1d*grid_points_1d; i++)
	{
		grid[i] = 0.0;
	}

	double mesh_width = 1.0/((double)(grid_points_1d-1));
	
	for (int i = 0; i < grid_points_1d; i++)
	{
		// x-boundaries
		grid[i] = eval_init_func(0.0, ((double)i)*mesh_width);
		grid[i + ((grid_points_1d)*(grid_points_1d-1))] = eval_init_func(1.0, ((double)i)*mesh_width);
		// y-boundaries
		grid[i*grid_points_1d] = eval_init_func(((double)i)*mesh_width, 0.0);
		grid[(i*grid_points_1d) + (grid_points_1d-1)] = eval_init_func(((double)i)*mesh_width, 1.0);
	}
}

/**
 * initializes the right hand side, we want to keep it simple and
 * solve the Laplace equation instead of Poisson (-> b=0)
 *
 * @param b the right hand side
 */
void init_b(double* b)
{
	int size_x = grid_points_x;
	if(rank > 0){
		size_x += 2;
	}
	
	// set all points to zero
	for (int i = 0; i < size_x*grid_points_y; i++)
	{
		b[i] = 0.0;
	}
}

/**
 * copies data from one grid to another
 *
 * @param dest destination grid
 * @param src source grid
 */
void g_copy(double* dest, double* src)
{
	int size_x = grid_points_x;
	if(rank > 0){
		size_x += 2;
	}
	
	for (int i = 0; i < size_x*grid_points_y; i++)
	{
		dest[i] = src[i];
	}
}

/**
 * calculates the dot product of the two grids (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 *
 * @param grid1 first grid
 * @param grid2 second grid
 */
 
double g_dot_product(double* grid1, double* grid2)
{
	double tmp = 0.0;

	#pragma omp parallel for
	for (int i = 1; i <= grid_points_x; i++)
	{
		for (int j = 1; j < grid_points_y-1; j++)
		{
			tmp += (grid1[(i*grid_points_y)+j] * grid2[(i*grid_points_y)+j]);
		}
	}
	
	return tmp;
}

/**
 * scales a grid by a given scalar (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 *
 * @param grid grid to be scaled
 * @param scalar scalar which is used to scale to grid
 */
void g_scale(double* grid, double scalar)
{
	#pragma omp parallel for
	for (int i = 1; i <= grid_points_x; i++)
	{
		for (int j = 1; j < grid_points_y-1; j++)
		{
			grid[(i*grid_points_y)+j] *= scalar;
		}
	}
}

/**
 * implements BLAS's Xaxpy operation for grids (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 *
 * @param dest destination grid
 * @param src source grid
 * @param scalar scalar to scale to source grid
 */
void g_scale_add(double* dest, double* src, double scalar)
{
	#pragma omp parallel for
	for (int i = 1; i <= grid_points_x; i++)
	{
		for (int j = 1; j < grid_points_y-1; j++)
		{
			dest[(i*grid_points_y)+j] += (scalar*src[(i*grid_points_y)+j]);
		}
	}
}

/**
 * implements the the 5-point finite differences stencil (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 * 
 * @param grid grid for which the stencil should be evaluated
 * @param result grid where the stencil's evaluation should be stored
 */
void g_product_operator(double* grid, double* result)
{
	
	// Send/Receive Border Columns
	if(rank > 1){
		MPI_Send(&grid[grid_points_y], 1, myVectorType1Col, rank-1, 0, MPI_COMM_WORLD);
		MPI_Recv(&grid[0], grid_points_y, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	if(rank < size-1){
		MPI_Send(&grid[grid_points_x*grid_points_y], 1, myVectorType1Col, rank+1, 0, MPI_COMM_WORLD);
		MPI_Recv(&grid[(1+grid_points_x)*grid_points_y], grid_points_y, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	
	double mesh_width = 1.0/((double)(grid_points_1d-1));

	#pragma omp parallel for
	for (int i = 1; i <= grid_points_x; i++)
	{
		for (int j = 1; j < grid_points_y-1; j++)
		{
			result[(i*grid_points_y)+j] =  (
							(4.0*grid[(i*grid_points_y)+j]) 
							- grid[((i+1)*grid_points_y)+j]
							- grid[((i-1)*grid_points_y)+j]
							- grid[(i*grid_points_y)+j+1]
							- grid[(i*grid_points_y)+j-1]
							) * (mesh_width*mesh_width);
		}
	}
}

/**
 * The CG Solver (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 *
 * For details please see :
 * http://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
 *
 * @param grid the grid containing the initial condition
 * @param b the right hand side
 * @param cg_max_iterations max. number of CG iterations 
 * @param cg_eps the CG's epsilon
 */
std::size_t solve(double* grid, double* b, std::size_t cg_max_iterations, double cg_eps)
{
	if(rank == 0){
		std::cout << "Starting Conjugated Gradients" << std::endl;
	}
		
	double eps_squared = cg_eps*cg_eps;
	std::size_t needed_iters = 0;
	
	if(rank == 0){
		// Separate Grid and send parts to all tasks
		int small_grid_points_x = (grid_points_1d-2)/(size-1);
		MPI_Type_vector(1, (small_grid_points_x+2)*grid_points_1d, 0, MPI_DOUBLE, &myVectorTypeGridWide);
		MPI_Type_commit(&myVectorTypeGridWide);
		for(int i=0;i<size-1;i++){
			MPI_Send(&grid[i*small_grid_points_x*grid_points_1d], 1, myVectorTypeGridWide, i+1, 0, MPI_COMM_WORLD);
		}
	}else{
		// Receive own Sub-Grid
		MPI_Recv(grid, (grid_points_x+2)*grid_points_y, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	double* q;
	double* r;
	double* d;
	double* b_save;
	
	if(rank > 0){
		// define temporal vectors
		q = (double*)_mm_malloc((grid_points_x+2)*grid_points_y*sizeof(double), 64);
		r = (double*)_mm_malloc((grid_points_x+2)*grid_points_y*sizeof(double), 64);
		d = (double*)_mm_malloc((grid_points_x+2)*grid_points_y*sizeof(double), 64);
		b_save = (double*)_mm_malloc((grid_points_x+2)*grid_points_y*sizeof(double), 64);
		
		g_copy(q, grid);
		g_copy(r, grid);
		g_copy(d, grid);
		g_copy(b_save, b);
	}
	
	double delta_0 = 0.0;
	double delta_old = 0.0;
	double delta_new = 0.0;
	double beta = 0.0;
	double a = 0.0;
	double residuum = 0.0;
	
	if(rank > 0){
		g_product_operator(grid, d);
		g_scale_add(b, d, -1.0);
		g_copy(r, b);
		g_copy(d, r);
	}

	// calculate starting norm
	// Do task, collect partial results and sum up
	if(rank > 0){
		delta_new = g_dot_product(r, r);
		MPI_Send(&delta_new, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}else{
		delta_new = 0;
		for(int i = 1; i < size; i++){
			double delta_buf = 0.0;
			MPI_Recv(&delta_buf, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			delta_new += delta_buf;
		}
	}
	// Send/Receive final result
	if(rank > 0){
		MPI_Recv(&delta_new, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}else{
		for(int i = 1; i < size; i++){
			MPI_Send(&delta_new, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
	}
	
	delta_0 = delta_new*eps_squared;
	residuum = (delta_0/eps_squared);
	
	if(rank == 0){
		std::cout << "Starting norm of residuum: " << (delta_0/eps_squared) << std::endl;
		std::cout << "Target norm:               " << (delta_0) << std::endl;
	}

	while ((needed_iters < cg_max_iterations) && (delta_new > delta_0))
	{
		// q = A*d
		if(rank > 0){
			g_product_operator(d, q);
		}

		// a = d_new / d.q
		// Do task, collect partial results and sum up
		if(rank > 0){
			double g_dot_product_var = g_dot_product(d, q);
			MPI_Send(&g_dot_product_var, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}else{
			double g_dot_product_tmp = 0.0;
			for(int i = 1; i < size; i++){
				double g_dot_product_buf = 0.0;
				MPI_Recv(&g_dot_product_buf, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				g_dot_product_tmp += g_dot_product_buf;
			}
			a = delta_new/g_dot_product_tmp;
		}
		// Send/Receive final result
		if(rank > 0){
			MPI_Recv(&a, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}else{
			for(int i = 1; i < size; i++){
				MPI_Send(&a, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}
		}
		
		if(rank > 0){
			// x = x + a*d
			g_scale_add(grid, d, a);
			
			if ((needed_iters % 50) == 0)
			{
				g_copy(b, b_save);
				g_product_operator(grid, q);
				g_scale_add(b, q, -1.0);
				g_copy(r, b);
			}
			else
			{
				// r = r - a*q
				g_scale_add(r, q, -a);
			}
		}
		
		// calculate new deltas and determine beta
		delta_old = delta_new;		
		
		// Do task, collect partial results and sum up
		if(rank > 0){
			delta_new = g_dot_product(r, r);
			MPI_Send(&delta_new, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}else{
			delta_new = 0;
			for(int i = 1; i < size; i++){
				double delta_buf = 0.0;
				MPI_Recv(&delta_buf, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				delta_new += delta_buf;
			}
		}
		// Send/Receive final result
		if(rank > 0){
			MPI_Recv(&delta_new, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}else{
			for(int i = 1; i < size; i++){
				MPI_Send(&delta_new, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}
		}
		beta = delta_new/delta_old;

		// adjust d
		if(rank > 0){
			g_scale(d, beta);
			g_scale_add(d, r, 1.0);
		}
		
		residuum = delta_new;
		needed_iters++;
		if(rank == 0){
			std::cout << "(iter: " << needed_iters << ")delta: " << delta_new << std::endl;
		}
	}
	
	if(rank == 0){
		std::cout << "Number of iterations: " << needed_iters << " (max. " << cg_max_iterations << ")" << std::endl;
		std::cout << "Final norm of residuum: " << delta_new << std::endl;
		
		// combine received grid-parts to whole grid
		int small_grid_points_x = (grid_points_1d-2)/(size-1);
		for(int i=0;i<size-1;i++){
			MPI_Recv(&grid[(1+i*small_grid_points_x)*grid_points_1d], small_grid_points_x*grid_points_y, MPI_DOUBLE, i+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}else{
	
		_mm_free(d);
		_mm_free(q);
		_mm_free(r);
		_mm_free(b_save);
	
		// Send out final grid-part
		MPI_Type_vector(1, grid_points_x*grid_points_y, 0, MPI_DOUBLE, &myVectorTypeGrid);
		MPI_Type_commit(&myVectorTypeGrid);
		MPI_Send(&grid[grid_points_y], 1, myVectorTypeGrid, 0, 0, MPI_COMM_WORLD);
	}
}

/**
 * main application
 *
 * @param argc number of cli arguments
 * @param argv values of cli arguments
 */
int main(int argc, char* argv[])
{
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	// check if all parameters are specified
	if (argc != 4)
	{
		if(rank == 0){
			std::cout << std::endl;
			std::cout << "meshwidth" << std::endl;
			std::cout << "cg_max_iterations" << std::endl;
			std::cout << "cg_eps" << std::endl;
			std::cout << std::endl;
			std::cout << "example:" << std::endl;
			std::cout << "./app 0.125 100 0.0001" << std::endl;
			std::cout << std::endl;
		}
		
		return -1;
	}
	
	// read cli arguments
	double mesh_width = atof(argv[1]);
	size_t cg_max_iterations = atoi(argv[2]);
	double cg_eps = atof(argv[3]);

	// calculate grid points per dimension
	grid_points_1d = (std::size_t)(1.0/mesh_width)+1;
	
	if(rank == 0){
		grid_points_x = grid_points_1d;
		grid_points_y = grid_points_1d;
	}else{
		// !!! (grid_points_1d-2) % (size-1) == 0 !!!
		grid_points_x = (grid_points_1d-2)/(size-1);
		grid_points_y = grid_points_1d;
		if(rank == 1){
			printf("grid_points_x is %d\n", grid_points_x);
			printf("grid_points_y is %d\n", grid_points_y);
		}
	}
	
	// initialize the grid and rights hand side
	double* grid, b;
	
	if(rank == 0){
		// Init Main-Grid for rank 0
		printf("grid_points_1d is %d\n", grid_points_1d);
		double* grid = (double*)_mm_malloc(grid_points_1d*grid_points_1d*sizeof(double), 64);
		double* b = (double*)_mm_malloc(grid_points_1d*grid_points_1d*sizeof(double), 64);
		init_grid(grid);
		store_grid(grid, "initial_condition.gnuplot");
		init_b(b);
		store_grid(b, "b.gnuplot");
		// solve Poisson equation using CG method
		timer_start();
		solve(grid, b, cg_max_iterations, cg_eps);
		double time = timer_stop();
		store_grid(grid, "solution.gnuplot");
		std::cout << std::endl << "Needed time: " << time << " s" << std::endl << std::endl;
		_mm_free(grid);
		_mm_free(b);
	}else{
		// All other ranks work on subset
		// no initialization, just allocation
		MPI_Type_vector(1, grid_points_1d, 0, MPI_DOUBLE, &myVectorType1Col);
		MPI_Type_commit(&myVectorType1Col);
		double* grid = (double*)_mm_malloc((grid_points_x+2)*grid_points_y*sizeof(double), 64);
		double* b = (double*)_mm_malloc((grid_points_x+2)*grid_points_y*sizeof(double), 64);
		init_b(b);
		solve(grid, b, cg_max_iterations, cg_eps);
		_mm_free(grid);
		_mm_free(b);
	}
	
	MPI_Finalize();
	
	return 0;
}


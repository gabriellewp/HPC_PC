/* row wise gaus algorithm
 * pattern for practical course
 * -------------------------
 * autor: Markus Brenk
 * date: 2002-09-25
 * =================================================== */

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <limits>
#include "Stopwatch.h"

/** print a 3x3 matrix */
void print_matrix(char* name, double matrix[3][3]);
/** print a 3d vector */
void print_vector(char* name, double vec[3]);
/**
 *  initialisation: generates the following LGS:
 *  ( 3  1  1)       (5r)                         (r)
 *  ( 1  4  1) * X = (6r)        => solution  X = (r)
 *  ( 1  1  5)       (7r)                         (r)
 */
void init(double r, double A[3][3], double b[3], double x[3]);
/** check if x_d == r */
bool check(double r, double x[3]);
/** performs gauss elimination */
void gauss_elimination(double A[3][3], double b[3], double x[3]);

int main() {
  double A[NRHS][3][3];
  double b[NRHS][3];
  double x[NRHS][3];
  
  for (int r = 0; r < NRHS; ++r) {
    init(static_cast<double>(r)/NRHS, A[r], b[r], x[r]);
  }

  Stopwatch stopwatch;
  stopwatch.start();

  for (int r = 0; r < NRHS; ++r) {
    gauss_elimination(A[r], b[r], x[r]);
  }

  double time = stopwatch.stop();
  printf("Time: %lf us\n", time * 1.0e6);
  
  bool correct = true;
  for (int r = 0; r < NRHS; r++) {
    correct = correct && check(static_cast<double>(r)/NRHS, x[r]);
  }
  if (!correct) {
    printf("Incorrect code.\n");
  }
  
  return 0;
}

void init(double r, double A[3][3], double b[3], double x[3]) {
	int i,j;
	int n = 3;

	for (j=0;j<n;j++)
	{
		b[j]=r * ((float)(2*n-2)+(float)(j+1));
		A[j][j]=(float)(n-1)+(float)(j+1);
		x[j]=0.;

		for (i=j+1;i<n;i++) {
			A[i][j]=1.;
		}
		for (i=0;i<j;i++) {
			A[i][j]=1.;
		}
	}
}

bool check(double r, double x[3]) {
  for (int d = 0; d < 3; ++d) {
    if (fabs(x[d] - r*1.0) > 10.0*std::numeric_limits<double>::epsilon()) {
      return false;
    }
  }
  return true;
}

void gauss_elimination(double A[3][3], double b[3], double x[3]) {
	int n = 3;
	int i,j,k;

	for (i = 0; i < n; i++) {

		for (j = i+1; j < n; j++) {
		  	A[i][j] = A[i][j] / A[i][i];
		}
		b[i] = b[i] / A[i][i];

		for (j = i+1; j < n; j++) {
			double factor = A[j][i];
			for (k = i; k < n; k++) {
				A[j][k] = A[j][k] - A[i][k] * factor;
			}
			b[j] = b[j] - factor * b[i];
		}
	}


	for (i = n-1; i >= 0; i--) {
		x[i] = b[i];
		for(j = i+1; j < n; j++) {
			x[i] -= A[i][j] * x[j];
		}
	}
}

void print_matrix(char* name, double matrix[3][3]) {
	int i, j;
	printf("Matrix %s: \n", name);
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			printf(" %f ", matrix[i][j]);
		}
		printf(" ;\n");
	}
}

void print_vector(char* name, double vec[3]) {
	int i;
	printf("vector %s: \n", name);
	for (i=0;i<3;i++)
	{
		printf(" %f \n",vec[i]);
	}
}


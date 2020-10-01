# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <math.h>
# include <omp.h>
const int N = 2048;
double R[N][N];
double Ri[N][N];

void print_matrix(int n, double M[n][n]);

void compute_inverse(int r0, int c0, int n);

double evaluate(int n);

int main (void)
{
	int n = N, np;	
	double sum = 0;
	double start, total_time, norm;
/*
 create a well-conditioned random R
*/
	for ( int i = 0; i < n; i ++)
	{
		for ( int j = i; j < n; j ++)
		{
			srand48(time(0));
			R[i][j] = drand48();
			sum += R[i][j];
		}
	}

	for ( int i = 0; i < n; i ++)
	{
		R[i][i] += sum;
	}
//	print_matrix(n, R);
/*
 compute the inverse of the upper triangular matrix recursively
*/
	start = omp_get_wtime();
	compute_inverse(0, 0, n);
	total_time = omp_get_wtime() - start; 
//	print_matrix(n, Ri);
/*
 calculate the error in computing inverse
*/
	norm = evaluate(n);
	#pragma omp parallel
	{
		np = omp_get_num_threads();
	}
	printf("n = %d, p = %d, error = %e, time (sec) = %8.4f\n", N, np, norm, total_time);
		
//	print_matrix(n, I);
	
	return 0;
}

void compute_inverse(int r0, int c0, int n)
{
	if (n == 1)
/*
 base case of recursion
*/ 
	{
		if (R[r0][c0] > 0)
			Ri[r0][c0] = 1/R[r0][c0];
	}
	else
	{
		int n1 = n/2;
/*
 compute the inverse of the two diagonal blocks
*/
		compute_inverse(r0, c0, n1);
		compute_inverse(r0 + n1, c0 + n1, n1);	
		double M[n1][n1];
		int i, j, k;
/*
 compute the off-diagonal block by matrix multiplication
*/
		# pragma omp parallel shared( M, R, Ri, r0, c0, n1) private (i, j, k)		
		{ 
		# pragma omp for	
		for ( i = 0; i < n1; i ++)
		{
			for ( j = 0; j < n1; j ++)
			{
				M[i][j] = 0.0;
				for ( k = 0; k < n1; k ++)
				{
					M[i][j] = M[i][j] - Ri[r0 + i][c0 + k] * R[r0 + k][c0 + n1 + j];
				}
			}
				
		}
		# pragma omp for
        for ( i = 0; i < n1; i ++)
		{
			for ( j = 0; j < n1; j ++)
			{
				Ri[r0 + i][c0 + n1 + j] = 0.0;
				for ( k = 0; k < n1; k ++)
				{
					Ri[r0 + i][c0 + n1 + j] = Ri[r0 + i][c0 + n1 + j] + M[i][k] * Ri[r0 + n1 + k][c0 + n1 + j];
				}
			}
		}
	  }
	}			
}

double evaluate(int n)
{
/*
 compute the product of R and Ri
 compare with identity matrix
*/
	double I[n][n];
	double error = 0.0;
	for (int i = 0; i < n; i ++)
	{
		for (int j = 0; j < n; j ++)
		{
			I[i][j] = 0.0;
			for (int k = 0; k < n; k ++)
			{
				I[i][j] += Ri[i][k] * R[k][j];
			}
			if (i == j)
				error += pow(I[i][j] - 1, 2);
			else
				error += pow(I[i][j],2);	
		}
	}
	return sqrt(error);
}

void print_matrix(int n, double M[n][n])
{
	for ( int i = 0; i < n; i ++)
	{
		for ( int j = 0; j < n; j ++)
		{
			printf ("%.1e ",M[i][j]);
		}
		printf ("\n");
	}	
}


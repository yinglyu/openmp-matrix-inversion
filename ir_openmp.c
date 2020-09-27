# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <omp.h>
double R[16][16];
double Ri[16][16];

void print_matrix(int n, double M[n][n]);

void compute_inverse(int r0, int c0, int n);

double evaluate(int n);

int main (void)
{
	int n = 16;	
	double sum = 0;
	for ( int i = 0; i < n; i ++)
	{
		for ( int j = i; j < n; j ++)
		{
			R[i][j] = drand48();
			sum += R[i][j];
		}
	}

	for ( int i = 0; i < n; i ++)
	{
		R[i][i] += sum;
	}
//	print_matrix(n, R);
	compute_inverse(0, 0, n);
//	print_matrix(n, Ri);
	double norm = evaluate(n);
	printf("Error in computing inverse: %e\n",norm);
		
//	print_matrix(n, I);
	
	return 0;
}

void compute_inverse(int r0, int c0, int n)
{
	if (n == 1)
	{
		if (R[r0][c0] > 0)
			Ri[r0][c0] = 1/R[r0][c0];
	}
	else
	{
		int n1 = n/2;
		compute_inverse(r0, c0, n1);
		compute_inverse(r0 + n1, c0 + n1, n1);	
		double M[n1][n1];
		for ( int i = 0; i < n1; i ++)
		{
			for (int j = 0; j < n1; j ++)
			{
				M[i][j] = 0.0;
				for ( int k = 0; k < n1; k ++)
				{
					M[i][j] = M[i][j] - Ri[r0 + i][c0 + k] * R[r0 + k][c0 + n1 + j];
				}
			}
				
		}
        for ( int i = 0; i < n1; i ++)
		{
			for (int j = 0; j < n1; j ++)
			{
				Ri[r0 + i][c0 + n1 + j] = 0.0;
				for ( int k = 0; k < n1; k ++)
				{
					Ri[r0 + i][c0 + n1 + j] = Ri[r0 + i][c0 + n1 + j] + M[i][k] * Ri[r0 + n1 + k][c0 + n1 + j];
				}
			}
		}
	}			
}

double evaluate(int n)
{
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


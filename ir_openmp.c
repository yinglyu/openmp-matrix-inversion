# include <stdlib.h>
# include <stdio.h>
# include <omp.h>

int main (void)
{
	double r[16][16];
	double ri[16][16];
	int n = 16;
	
	double sum = 0;
	for ( int i = 0; i < n; i ++)
	{
		for ( int j = i; j < n; j ++)
		{
			r[i][j] = drand48();
			sum += r[i][j];
		}
	}

	for ( int i = 0; i < n; i ++)
	{
		r[i][i] += sum;
	}

	for ( int i = 0; i < n; i ++)
	{
		for ( int j = 0; j < n; j ++)
		{
			printf ("%.2f ",r[i][j]);
		}
		printf ("\n");
	}
	return 0;
}


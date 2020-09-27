# include <stdlib.h>
# include <stdio.h>
# include <omp.h>

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


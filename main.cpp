#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include "TJunction.h"
#include <omp.h>

// Calculates the S-parameters for a T-Junction waveguide for 21
// different widths of the connecting aperture and prints them to
// the terminal. Accepts two optional arguments: the number of 
// basis functions and the number of waveguide modes to include

/*int main()
{
	TJunction junction(5e9, 47.55e-3, 22.15e-3, 10e-3, 5, 100);
	double* S = junction.SolveTJunction();
	printf("S11: %f\nS22: %f\nS33: %f\n\n", S[0], S[1], S[2]);

	return 0;
}*/

int main(int argc, char *argv[])
{	
	double a = 47.55e-3;
	double b = 22.15e-3;
	double f = 5e9;

	long int basisFunctions, maxModeNumber;

	if (argc == 1) // No input arguments, set to default
	{
		basisFunctions = 5;
		maxModeNumber = 100;
	}
	else if (argc == 3)
	{
		basisFunctions = atoi(argv[1]);
		maxModeNumber = atoi(argv[2]);
	}
	else // Wrong number of arguments!
	{
		printf("\nIllegal number of arguments! Either zero or two, where the inputs correspond to the number of basis functions and waveguide modes, respectively\n\n");
		return 0;
	}

	printf("\nw/a \t\t S11 \t\t S21 \t\t S31\n\n");

	TJunction waveguide;
	double start_time = (double)clock();
	#pragma omp parallel for
	for (int i = 0; i <= 20; i++)
	{
		// Calculate S-parameters and print them
		double* SParameters = waveguide.SolveTJunction(f, a, b, a/20.0*i, basisFunctions, maxModeNumber);
		printf("%f \t %f \t %f \t %f \n", 1.0/20.0*i, SParameters[0], SParameters[1], SParameters[2]);
	}
	double end_time = (double)clock();
	double time_lapsed = (end_time - start_time) / (double)CLOCKS_PER_SEC;

	printf("\nTime elapsed: %f s\n\n", time_lapsed);

	return 0;
}


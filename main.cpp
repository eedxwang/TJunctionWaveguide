#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include "SolveTJunction.h"

// Calculates the S-parameters for a T-Junction waveguide for 21
// different widths of the connecting aperture and prints them to
// the terminal. Accepts two optional arguments: the number of 
// basis functions and the number of waveguide modes to include
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

	double start_time = (double)clock();
	for (double w = 0.0; w <= a; w += 0.05*a)
	{
		double* SParameters = SolveTJunction(f, a, b, w, basisFunctions, maxModeNumber); // Calculate S-parameters
		printf("%f \t %f \t %f \t %f \n", w/a, SParameters[0], SParameters[1], SParameters[2]);
	}
	double end_time = (double)clock();
	double time_lapsed = (end_time - start_time) / (double)CLOCKS_PER_SEC;

	printf("\nTime elapsed: %f s\n\n", time_lapsed);

	return 0;
}


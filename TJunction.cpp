#include <complex>
#include "f2c.h"
#include "clapack.h"
#include <stdio.h>
#include "TJunction.h"

#define PI 3.14159265358979323846
#define ETA 376.73			// Impedance of free space

// Returns an array of doubles, containing the magnitude of the scattering parameters S11, S21 and S31 in dB
double* TJunction::SolveTJunction(double freq, double a, double b, double w, long int basisFunctions, long int maxModeNumber)
{
	if (w == 0)				// This is illegal due to the occurence of w in denominators
		w = 0.00001*a;

	integer N = basisFunctions;

	double c = 299792458.0;				// Speed of light
	double k_0 = 2.0 * PI*freq / c;		// Free space wavenumber 

	double d = 3.0 / 4.0 * (2.0 * PI) / (std::sqrt(k_0*k_0 - PI*PI / (a*a))); // Length of virtual cavity (= 3/4*lambda_prop)

	// Create lookup tables
	std::complex<double>* gamma_u_a = NULL;
	std::complex<double>* gamma_u_d = NULL;
	std::complex<double>* factor1_lookup = NULL;
	std::complex<double>* iZplusTable = NULL;
	std::complex<double>* iZminusTable = NULL;
	std::complex<double>* cothLookup = NULL;
	createGammaULookupTables(a, d, k_0, maxModeNumber, &gamma_u_a, &gamma_u_d);
	createFactor1LookupTable(a, b, d, k_0, maxModeNumber, gamma_u_a, &factor1_lookup);
	createIzLookup(a, w, N, maxModeNumber, gamma_u_a, &iZplusTable, &iZminusTable);
	createCothLookup(a, b, d, k_0, maxModeNumber, gamma_u_d, &cothLookup);

	// Create matrices for LAPACK
	doublecomplex* A = (doublecomplex*)malloc(N*N*sizeof(doublecomplex));
	doublecomplex* RHS = (doublecomplex*)malloc(N*sizeof(doublecomplex));
	std::complex<double> matrixElement;

	for (int i = 0; i < N; i++)
	{
		matrixElement = rightHandSide(a, b, w, k_0, i + 1, maxModeNumber, gamma_u_a, iZminusTable);
		*(RHS + i) = { std::real(matrixElement), std::imag(matrixElement) };
		for (int j = i; j < N; j++)
		{
			matrixElement = selfAdmittanceRegion1And2(a, b, w, k_0, d, i + 1, j + 1, maxModeNumber,
				gamma_u_a, gamma_u_d, factor1_lookup, iZplusTable, iZminusTable, cothLookup);
			*(A + N*j + i) = { std::real(matrixElement), std::imag(matrixElement) };
			*(A + N*i + j) = *(A + N*j + i);						// Due to reciprocity/symmetry of the matrix
		}
	}

	// Solve system
	integer* ipiv = (integer*)malloc(N*sizeof(integer));
	integer INFO = 0;
	integer numRHS = 1;

	zgesv_(&N, &numRHS, A, &N, ipiv, RHS, &N, &INFO);

	// Check if a unique solution to the system was found
	if (INFO != 0)
	{
		printf("No unique solution to the system was found!\n INFO=%f\n Exiting... \n\n", INFO);
		free(A);
		free(RHS);
		return 0;
	}

	// Extract the expansion coefficients from the solution and put them in the (complex) vector M
	std::complex<double>* M = (std::complex<double>*)malloc(N*sizeof(std::complex<double>));
	for (int n = 0; n < N; n++)
	{
		*(M + n) = std::complex<double>((*(RHS + n)).r, (*(RHS + n)).i);
	}

	// Calculate scattering parameters
	std::complex<double> S11 = calculateS11(a, b, w, N, k_0, M, maxModeNumber, gamma_u_a, iZminusTable);
	std::complex<double> S21 = calculateS21(a, b, w, N, k_0, M, maxModeNumber, gamma_u_a, iZplusTable);
	std::complex<double> S31 = calculateS31(a, b, w, N, k_0, M, gamma_u_a);

	// Convert magnitudes to dB
	double absS11 = (std::real(S11))*(std::real(S11)) + (std::imag(S11))*(std::imag(S11));
	double S11_dB = 10 * std::log10(absS11);
	double absS21 = (std::real(S21))*(std::real(S21)) + (std::imag(S21))*(std::imag(S21));
	double S21_dB = 10 * std::log10(absS21);
	double absS31 = (std::real(S31))*(std::real(S31)) + (std::imag(S31))*(std::imag(S31));
	double S31_dB = 10 * std::log10(absS31);

	// Return dB values
	double* scatteringParameters = (double*)malloc(3 * sizeof(double));
	scatteringParameters[0] = S11_dB;
	scatteringParameters[1] = S21_dB;
	scatteringParameters[2] = S31_dB;

	free(A);
	free(RHS);
	free(gamma_u_a);
	free(gamma_u_d);
	free(factor1_lookup);
	free(iZplusTable);
	free(iZminusTable);
	free(cothLookup);

	return scatteringParameters;
}

// Calculates the integral A.2 in the thesis. SEEMS TO WORK CORRECTLY
double TJunction::integralXCoupling(int n, int k, double a, double w)
{
	// First check if length is legal:
	if (abs(n*n*PI*PI / (a*a) - k*k*PI*PI / (w*w)) < 1e-10)
	{
		// Illegal length! Change w sliiiightly
		w = w + 0.00001*a;
	}

	double result = k*PI / w / (n*n*PI*PI / (a*a) - k*k*PI*PI / (w*w))*(std::pow(-1, k)*std::sin(n*PI / 2 + n*PI*w / (2.0*a)) - std::sin(n*PI / 2.0 - n*PI*w / (2.0*a)));
	return result;
}

// Calculates the integral A.1 in the thesis SEEMS TO WORK CORRECTLY
std::complex<double> TJunction::integralZCoupling(int n, int k, double a, double w, std::complex<double> gamma_u, int sign)
{
	std::complex<double> result = k*PI / w*std::exp(-0.5*sign*gamma_u*w) / (std::pow(gamma_u, 2) + k*k*PI*PI / (w*w))*(1.0 - (((double)std::pow(-1, k))*std::exp((1.0*sign)*gamma_u*w)));
	return result;
}

// SEEMS TO WORK CORRECTLY!
std::complex<double> TJunction::selfAdmittanceRegion1And2(double a, double b, double w, double k_0, double d, int i, int j, int maxModeNumber,
	std::complex<double>* gamma_u_a, std::complex<double>* gamma_u_d, std::complex<double>* factor1Lookup,
	std::complex<double>* iZplusTable, std::complex<double>* iZminusTable, std::complex<double>* cothLookup)
{
	std::complex<double> selfAdmittance(0.0, 0.0);
	std::complex<double> AdmIX(0.0, 0.0);
	std::complex<double> AdmIZ(0.0, 0.0);
	std::complex<double> AdmII(0.0, 0.0);
	std::complex<double> gamma_u;						// Propagation constant
	std::complex<double> root;
	std::complex<double> I(0.0, 1.0);					// Imaginary unit
	std::complex<double> newTerm, factor1, factor2;

	for (int mode = 1; mode <= maxModeNumber; mode++)
	{
		// Calculate contribution from x-travelling modes
		gamma_u = gamma_u_d[mode - 1];
		newTerm = -cothLookup[mode - 1] * integralXCoupling(mode, i, d, w)*integralXCoupling(mode, j, d, w);
		AdmIX += newTerm;

		// Calculate contribution from z-travelling modes
		gamma_u = gamma_u_a[mode - 1];
		factor1 = factor1Lookup[mode - 1];
		factor2 = iZminusTable[mode - 1 + (i - 1)*maxModeNumber] * iZminusTable[mode - 1 + (j - 1)*maxModeNumber] -
			std::exp(-gamma_u*d)*iZplusTable[mode - 1 + (i - 1)*maxModeNumber] * iZminusTable[mode - 1 + (j - 1)*maxModeNumber] +
			iZplusTable[mode - 1 + (i - 1)*maxModeNumber] * iZplusTable[mode - 1 + (j - 1)*maxModeNumber] -
			std::exp(-gamma_u*d)*iZminusTable[mode - 1 + (i - 1)*maxModeNumber] * iZplusTable[mode - 1 + (j - 1)*maxModeNumber];
		newTerm = factor1*factor2;
		AdmIZ += newTerm;

		// Calculate contribution in region II
		newTerm = -2.0 * b*gamma_u / (a*I*k_0*ETA)*integralXCoupling(mode, i, a, w)*integralXCoupling(mode, j, a, w);
		AdmII += newTerm;
	}
	selfAdmittance = AdmIX + AdmIZ + AdmII;

	return selfAdmittance;
}

// Calculates element i in the right hand side vector
std::complex<double> TJunction::rightHandSide(double a, double b, double w, double k_0, int i, long int maxModeNumber,
	std::complex<double>* gamma_u_a, std::complex<double>* iZminusLookup)
{
	double k_cut = PI / a;
	std::complex<double> I(0.0, 1.0);
	std::complex<double> gamma_u = gamma_u_a[0];

	std::complex<double> RHS = k_cut / gamma_u*std::sqrt(2.0 * gamma_u*b / (I*k_0*ETA*a))*iZminusLookup[(i - 1)*maxModeNumber];

	return RHS;
}

void TJunction::printMatrix(integer rows, integer cols, doublecomplex* matrix)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			printf("%e+%ei\t", (*(matrix + rows*j + i)).r, (*(matrix + rows*j + i)).i);
		}
		printf("\n");
	}
}

std::complex<double> TJunction::calculateS11(double a, double b, double w, long int N, double k_0, std::complex<double>* M, long int maxModeNumber,
	std::complex<double>* gamma_u_a, std::complex<double>* iZminusLookup)
{
	std::complex<double> S11(0.0, 0.0);
	std::complex<double> I(0.0, 1.0);
	std::complex<double> gamma_u = gamma_u_a[0];

	for (int i = 0; i < N; i++)
	{
		S11 += (*(M + i))*iZminusLookup[i*maxModeNumber];
	}
	S11 = S11*0.5*PI / (a*gamma_u)*std::sqrt(2.0 * b*gamma_u / (a*I*k_0*ETA));

	return S11;
}

std::complex<double> TJunction::calculateS21(double a, double b, double w, long int N, double k_0, std::complex<double>* M,
	long int maxModeNumber, std::complex<double>* gamma_u_a, std::complex<double>* iZplusLookup)
{
	std::complex<double> S21(0.0, 0.0);
	std::complex<double> I(0.0, 1.0);
	std::complex<double> gamma_u = gamma_u_a[0];

	for (int i = 0; i < N; i++)
	{
		S21 += (*(M + i))*iZplusLookup[maxModeNumber*i];
	}
	S21 = S21*PI / (2.0*a)*sqrt(2.0*b / (gamma_u*I*k_0*ETA*a)) + 1.0;

	return S21;
}

std::complex<double> TJunction::calculateS31(double a, double b, double w, long int N, double k_0, std::complex<double>* M, std::complex<double>* gamma_u_a)
{
	std::complex<double> S31(0.0, 0.0);
	std::complex<double> I(0.0, 1.0);
	std::complex<double> gamma_u = gamma_u_a[0];

	for (int i = 0; i < N; i++)
	{
		S31 += (*(M + i))*integralXCoupling(1, i + 1, a, w);
	}
	S31 = S31*(-1.0*sqrt(2.0 * gamma_u*b / (a*I*k_0*ETA)));

	return S31;
}

// Creates lookup table for all gamma_u from n=1 to n=maxModeNumber for lengths a and d
void TJunction::createGammaULookupTables(double a, double d, double k_0, long int maxModeNumber,
	std::complex<double>** gamma_u_a, std::complex<double>** gamma_u_d)
{
	std::complex<double>* gamma_u_a_table = (std::complex<double>*)malloc(maxModeNumber*sizeof(std::complex<double>));
	std::complex<double>* gamma_u_d_table = (std::complex<double>*)malloc(maxModeNumber*sizeof(std::complex<double>));

	std::complex<double> root;
	for (int n = 1; n <= maxModeNumber; n++)
	{
		root = std::complex<double>(PI*PI*n*n / (a*a) - k_0*k_0, 0.0);
		gamma_u_a_table[n - 1] = std::sqrt(root);

		root = std::complex<double>(PI*PI*n*n / (d*d) - k_0*k_0, 0.0);
		gamma_u_d_table[n - 1] = std::sqrt(root);
	}

	*gamma_u_a = gamma_u_a_table;
	*gamma_u_d = gamma_u_d_table;
}

// Create lookup table useful for calculating contribution from z-modes to self admittance
// in region I
void TJunction::createFactor1LookupTable(double a, double b, double d, double k_0, long int maxModeNumber, std::complex<double>* gamma_u_a,
	std::complex<double>** factor1Lookup)
{
	std::complex<double> gamma_u;
	std::complex<double> I(0.0, 1.0);
	std::complex<double>* lookup_table = (std::complex<double>*)malloc(maxModeNumber*sizeof(std::complex<double>));

	for (int n = 1; n <= maxModeNumber; n++)
	{
		gamma_u = gamma_u_a[n - 1];
		lookup_table[n - 1] = n*n*PI*PI*b / (2 * a*a*a*gamma_u*I*k_0*ETA*std::sinh(gamma_u*d));
	}

	*factor1Lookup = lookup_table;
}

// Create lookup tables for z-integrals
void TJunction::createIzLookup(double a, double w, long int N, long int maxModeNumber, std::complex<double>* gamma_u_a,
	std::complex<double>** iZplusLookup, std::complex<double>** iZminusLookup)
{
	std::complex<double>* iz_plus_table = (std::complex<double>*)malloc(N*maxModeNumber*sizeof(std::complex<double>));
	std::complex<double>* iz_minus_table = (std::complex<double>*)malloc(N*maxModeNumber*sizeof(std::complex<double>));

	int nIndex, kIndex;
	for (int n = 1; n <= maxModeNumber; n++)
	{
		for (int k = 1; k <= N; k++)
		{
			nIndex = n - 1;
			kIndex = k - 1;
			iz_plus_table[nIndex + kIndex*maxModeNumber] = integralZCoupling(n, k, a, w, gamma_u_a[nIndex], 1);
			iz_minus_table[nIndex + kIndex*maxModeNumber] = integralZCoupling(n, k, a, w, gamma_u_a[nIndex], -1);
		}
	}

	*iZplusLookup = iz_plus_table;
	*iZminusLookup = iz_minus_table;
}

void TJunction::createCothLookup(double a, double b, double d, double k_0, long int maxModeNumber, std::complex<double>* gamma_u_d,
	std::complex<double>** cothLookup)
{
	std::complex<double>* lookup = (std::complex<double>*)malloc(maxModeNumber*sizeof(std::complex<double>));
	std::complex<double> I(0.0, 1.0);
	std::complex<double> gamma_u;

	for (int n = 1; n <= maxModeNumber; n++)
	{
		gamma_u = gamma_u_d[n - 1];
		lookup[n - 1] = std::cosh(gamma_u*a) / std::sinh(gamma_u*a)*2.0*gamma_u*b / (I*k_0*ETA*d);
	}

	*cothLookup = lookup;
}
#pragma once

#include <complex>
#include "f2c.h"
#include "clapack.h"
#include <stdio.h>

// This is the main function that should be called to calculate the S-parameters for a T-Junction waveguide,
// the rest of the functions are simply support functions that are called by this function
double* SolveTJunction(double freq, double a, double b, double w, long int basisFunctions, long int maxModeNumber);

// Other functions
double integralXCoupling(int n, int k, double a, double w);

std::complex<double> integralZCoupling(int n, int k, double a, double w, std::complex<double> gamma_u, int sign);

std::complex<double> selfAdmittanceRegion1And2(double a, double b, double w, double k_0, double d, int i, int j, int maxModeNumber,
	std::complex<double>* gamma_u_a, std::complex<double>* gamma_u_d, std::complex<double>* factor1Lookup,
	std::complex<double>* iZplusTable, std::complex<double>* iZminusTable, std::complex<double>* cothLookup);

std::complex<double> rightHandSide(double a, double b, double w, double k_0, int i, long int maxModeNumber,
	std::complex<double>* gamma_u_a, std::complex<double>* iZminusLookup);

void printMatrix(integer rows, integer cols, doublecomplex* matrix);

std::complex<double> calculateS11(double a, double b, double w, long int N, double k_0, std::complex<double>* M, long int maxModeNumber,
	std::complex<double>* gamma_u_a, std::complex<double>* iZminusLookup);

std::complex<double> calculateS21(double a, double b, double w, long int N, double k_0, std::complex<double>* M,
	long int maxModeNumber, std::complex<double>* gamma_u_a, std::complex<double>* iZplusLookup);

std::complex<double> calculateS31(double a, double b, double w, long int N, double k_0, std::complex<double>* M, std::complex<double>* gamma_u_a);

void createGammaULookupTables(double a, double d, double k_0, long int maxModeNumber,
	std::complex<double>** gamma_u_a, std::complex<double>** gamma_u_d);

void createFactor1LookupTable(double a, double b, double d, double k_0, long int maxModeNumber, std::complex<double>* gamma_u_a,
	std::complex<double>** factor1Lookup);

void createIzLookup(double a, double w, long int N, long int maxModeNumber, std::complex<double>* gamma_u_a,
	std::complex<double>** iZplusLookup, std::complex<double>** iZminusLookup);

void createCothLookup(double a, double b, double d, double k_0, long int maxModeNumber, std::complex<double>* gamma_u_d,
	std::complex<double>** cothLookup);

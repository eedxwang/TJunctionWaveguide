The TJunctionProgram calculates the S-parameters for a T-Junction waveguide using a moment method
approach based on Seki's alternate expression for the waveguide Green's function. LAPACK (or,
specifically, CLAPACK) is used to solve the resulting matrix equation. 

The S-parameters are calculated and printed for 20 different widths of the connecting aperture.

The program can take two optional input arguments: the number of basis functions to use
in the expansion of the equivalent magnetic current, and the number of waveguide modes to
include in the modal expansion for the Green's function. 

Example call: 
>> TJunctionWaveguide.exe 5 100

It seems the agreement with FEKO is quite good.
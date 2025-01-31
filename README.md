# ExplicitlyCorrelatedGaussians_StochasticGrowth
C++ for stochastic optimization of wavefunctions for 2-4 electron atoms in a basis set of explicitly correlated gaussian functions. 
- Minimizes the energy for 2-4 electron systems in a basis set of spherical explicitly correlated gaussian functions [see, for example, 5.	M. Cafiero and L. Adamowicz, "Analytical Gradients for Singer's Molecular n-electron Explicitly Correlated Gaussians", Int. J. Quantum. Chem., 82, 151, (2001)].
- Stochastic optimization of the basis set perfomed by randomly modifying basis function coefficients of any basis function with a c^2 coefficient > 0.9, and then keeping the result only if the energy decreased.
- Checks for linear dependence by finding small eigenvalues of the overlap matrix.
- Can read in a basis set guess (necg.cpp) or generate one randomly (necg_grow.cpp) or read in a partial basis set and grow the rest (necg_read_grow.cpp). 
- Input: see attached examples. Must contain nuclear charge, number of electrons, and an initial basis set guess (necg.cpp), or nuclear charge, number of electrons, maximum number of basis functions to generate, and maximum sub-iterations per basis function (necg_grow.cpp), or nuclear charge, number of electrons, maximum number of basis functions to generate, maximum sub-iterations per basis function, and an initial basis set guess (necg_read_grow.cpp).
- Input files must be called "coeffs_in.csv." Examples for grow and read_grow must be changed to this name before use. 
- Output: A file with information on possible linear dependencies (generated each time one is found), a file with all eigenvalues sorted from lowest to highest (generated every iteration), and a file with the basis set coefficients (generated every 10 iterations, and at the end).
- Also includes code to genrate a random basis set outside of the main code (gen_basis_set.cpp). This takes no input file: all information is taken at the command line.
- Compiling: C++ 14, uses the Eigen library (https://eigen.tuxfamily.org/index.php?title=Main_Page). This library mustbe included in the compiler's path. 

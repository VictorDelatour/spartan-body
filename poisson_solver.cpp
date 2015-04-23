//#include <iostream>
#include <stdio.h>
#include "gravity_solver.hpp"

///
/// @file poisson_solver.cpp
/// @brief Fonction accessed from Fortran to solve the Poisson equation for the gravitiational potential
///     @author Hugo Babel - EPFL.
///     @date 2015-04-17 (initial version)


extern "C" void poisson_solver(int* nx, int* ny, int* nz, double* rho);

void poisson_solver(int* nx, int* ny, int* nz, double* rho)
{
    //std::cout << "\n This is a C++ code\n" << std::endl;
	printf("\nStarting to work in C++ \n");
	printf("Setup for solving poisson \n");
	printf("Grid has size %d x %d x %d \n\n", *nx, *ny, *nz);
	
	int temp(2);

	for(int i(0); i < temp; ++i){
		for(int j(0); j < temp; ++j){
			printf("rho[%d, %d, 1] = %f \n", i, j, rho[i + j*(*nx)]);
		}
	}
	printf("\n");
	
	GravitySolver solver(nx, ny, nz, rho);
	
	solver.init();
	
	solver.set_multiwavelet_order(3);
	solver.set_threshold(1e-1);
	
	printf("Projecting density\n");
	solver.project_density();
	printf("Density projected\n\n");
	
	
	printf("Printing density\n");
	solver.print_density(128);
	printf("Density printed\n\n");
	
	// printf("Solving potential\n");
	// solver.solve_potential(1e-2, 1e-4);
	// printf("Potential solved\n\n");
	//
	// printf("Printing potential\n");
	// solver.print_potential(128);
	// printf("Potential printed\n\n");

	
	
}
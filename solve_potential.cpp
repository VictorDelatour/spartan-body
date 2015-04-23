#include "solve_potential.hpp"

void solve_potential(World *world, const int& nx, const int& ny, const int& nz, real_t* density){
	
	
	real_function_3d rho_interp;
	real_function_3d phi;
	
	if (world->rank() == 0) printf("Setup initial parameters\n");
	set_initial_parameters(nx);
	if (world->rank() == 0) printf("Set...\n\n");
	
	if (world->rank() == 0) printf("Setup projection precision\n");
	set_projection_precision(4, 1e-2);
	if (world->rank() == 0) printf("Set...\n\n");
	
	if (world->rank() == 0) printf("Build projected density\n");
	build_projected_density(world, nx, ny, nz, density, rho_interp);
	if (world->rank() == 0) printf("Built...\n\n");
	
	if (world->rank() == 0) printf("Printing\n");
	// build_projected_density(world, nx, ny, nz, density, rho_interp);
	print_density(world, rho_interp, 128, nx);
	if (world->rank() == 0) printf("Printed...\n\n");
}

void set_initial_parameters(const int& nx){
	
	BoundaryConditions<3> bc(BC_PERIODIC);
	
	FunctionDefaults<3>::set_cubic_cell((double) 1, (double) nx);
	FunctionDefaults<3>::set_bc(bc); 
	FunctionDefaults<3>::set_apply_randomize(true);
	FunctionDefaults<3>::set_autorefine(true);
	FunctionDefaults<3>::set_refine(true);
	
}

void set_projection_precision(const int& order, const double& threshold){
	
	FunctionDefaults<3>::set_k(order);
	FunctionDefaults<3>::set_thresh(threshold);
	
}

void build_projected_density(World *world, const int& nx, const int& ny, const int& nz, real_t* density, real_function_3d& projected_density){
	
	real_functor_3d density_functor;
	
	density_functor = real_functor_3d(new DensityProjector(nx, ny, nz, &density[0]));
	projected_density = real_factory_3d(*world).functor(density_functor);
	
}

void print_density(World *world, const real_function_3d& projected_density, const int& numpts, const int& nx){
	
	// Check if the density has been projected!
	
	const char filename_density[] = "data/spartan_density.vts";
	
	Vector<double, 3> plotlo, plothi;
	Vector<long, 3> npoints;
	
	
	for(int i(0); i < 3; ++i){
		plotlo[i] = 1;
		plothi[i] = (double) nx;
		npoints[i] = numpts;
	}
	
	plotvtk_begin(*world, filename_density, plotlo, plothi, npoints);
	plotvtk_data(projected_density, "density", *world, filename_density, plotlo, plothi, npoints);
	plotvtk_end<3>(*world, filename_density);
	
}

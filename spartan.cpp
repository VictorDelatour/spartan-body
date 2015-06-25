#include <vector>
#include <cstdio>
#include <cstring>
#include <chrono>
#include <iostream>
#include <fstream>

#include <madness/mra/mra.h>

#include "density_projector3d.hpp"

using namespace madness;

typedef double real_t;
typedef Vector<double,3> coordT;

// Fortran functions subroutines called by spartan

extern "C" void get_dim_(int* nx, int* ny, int* nz, int* nparticles, int* nproc);
extern "C" void part_init_(const int* nx, const int* ny, const int* nz, const int* nparticles, const int* nproc, real_t* x, real_t* y, real_t* z, real_t* vx, real_t* vy, real_t* vz, real_t* mass);
extern "C" void project_density_(const int* nx, const int* ny, const int* nz, const int* nparticles, const real_t* x, const real_t* y, const real_t* z, real_t* mass, real_t* density, const int* step);

const double p( 1/( 4 * constants::pi ) );//*6.67384e-11));

///
/// @brief Initialization of madness parameters
///
/// @param[in]	nx	Upper limit of the box
///
/// Initializes madness parameters, with periodic BC, cubic cell and refinement of the mesh
///

void set_initial_parameters(const int& nx){
	
	BoundaryConditions<3> bc(BC_PERIODIC);
	
	FunctionDefaults<3>::set_cubic_cell(1.0, static_cast<double>(nx));
	FunctionDefaults<3>::set_bc(bc);
	FunctionDefaults<3>::set_apply_randomize(true);
	FunctionDefaults<3>::set_autorefine(true);
	FunctionDefaults<3>::set_refine(true);

}

///
/// @brief Set precision for the multiwavelet representation
///
/// @param[in]	order		Order of the wavelets used
/// @param[in]	threshold	Accuracy
///


void set_projection_precision(const int& order, const double& threshold){
	
	FunctionDefaults<3>::set_k(order);
	FunctionDefaults<3>::set_thresh(threshold);
	
}

///
/// @brief Print density to vtk file
///
/// @param[in]	world				Reference to the madness:world
/// @param[in]	projected_density	Density to be printed
/// @param[in]	numpts				Number of points for the plot
/// @param[in]	nx					Upper limit of the plot box
///


void print_density(World& world, const real_function_3d& projected_density, const int& numpts, const int& nx){
	
	const char filename_density[] = "data/spartan_density.vts";
	
	Vector<double, 3> plotlo, plothi;
	Vector<long, 3> npoints;
	
	for(int i(0); i < 3; ++i){
		plotlo[i] = 1;
		plothi[i] = (double) nx;
		npoints[i] = numpts;
	}
	
	plotvtk_begin(world, filename_density, plotlo, plothi, npoints);
	plotvtk_data(projected_density, "density", world, filename_density, plotlo, plothi, npoints);
	plotvtk_end<3>(world, filename_density);
	
}

///
/// @brief Print potential to vtk file
///
/// @param[in]	world		Reference to the madness:world
/// @param[in]	potential	Potential to be printed
/// @param[in]	numpts		Number of points for the plot
/// @param[in]	nx			Upper limit of the plot box
///

void print_potential(World& world, const real_function_3d& potential, const int& numpts, const int& nx){
	
	const char filename_potential[] = "data/spartan_potential.vts";
	
	Vector<double, 3> plotlo, plothi;
	Vector<long, 3> npoints;
	
	
	for(int i(0); i < 3; ++i){
		plotlo[i] = 1;
		plothi[i] = (double) nx;
		npoints[i] = numpts;
	}
	
	plotvtk_begin(world, filename_potential, plotlo, plothi, npoints);
	plotvtk_data(potential, "potential", world, filename_potential, plotlo, plothi, npoints);
	plotvtk_end<3>(world, filename_potential);
	
}

///
/// @brief Built multiwavelet representation of the density
///
/// @param[in]		world				Reference to the madness:world
/// @param[in]		nx					x-dimension of the mesh
/// @param[in]		ny					y-dimension of the mesh
/// @param[in]		nz					z-dimension of the mesh
/// @param[in, out]	density				Array (size nx*ny*nz) containing the projection of the particles on the uniform
/// @param[in, out]	projected_density	Function to be initialized
///


void build_projected_density(World& world, const int& nx, const int& ny, const int& nz, real_t* density, real_function_3d& projected_density){
	
	if (world.rank() == 0) printf("\nDensity step\n");
	
	// Compute the coefficients for the cubic B-spline interpolation
	real_functor_3d density_functor = real_functor_3d(new DensityProjector(nx, ny, nz, &density[0]));
	

	// Do the cubic B-spline interpolation
	projected_density = real_factory_3d(world).functor(density_functor);

	
	print_density(world, projected_density, 128, 128);
	
}

///
/// @brief Compute gravitational potential from density
///
/// @param[in]		world				Reference to the madness:world
/// @param[in]		projected_density	Multiwavelet representation of the density
/// @param[in, out]	potential			S
/// @param[in]		precision			Accuracy of the coulomb operator
/// @param[in]		threshold			Accuracy of the coulomb operator
///
///

void compute_potential(World& world, const real_function_3d& projected_density, real_function_3d& potential, const double& precision, const double& threshold){
	
	if (world.rank() == 0) printf("\nPotential step\n");
	

	
	real_convolution_3d coulomb_operator = CoulombOperator(world, precision, threshold);
	
	potential = coulomb_operator(projected_density);
		
	// // Set mean of potential to 0
	// double integral = potential.norm2();
	// double volume = FunctionDefaults<3>::get_cell_volume();
	// double mean = integral*integral/volume;
	// potential = potential - mean;

}

///
/// @brief Solve Poisson equation
///
/// @param[in]		world			Reference to the madness:world
/// @param[in]		x				x-coordinate of the particles
/// @param[in]		y				y-coordinate of the particles
/// @param[in]		z				z-coordinate of the particles
/// @param[in]		nx				x-dimension of the cubic uniform mesh
/// @param[in]		ny				y-dimension of the cubic uniform mesh
/// @param[in]		nz				z-dimension of the cubic uniform mesh
/// @param[in]		nparticles		Number of particles used
/// @param[in, out]	density			Cubic uniform mesh containing the projection of the density from the particles
///
/// Main function of the code, solves the Poisson equation. It starts by setting initial parameters (BC, accuracy),
/// build a multiwavelet representation of the density and then solves the Poisson equation for the gravitational
/// potential using the projection of the density.
///

real_function_3d solve_potential(World& world, real_t* x, real_t* y, real_t* z, const int& nx, const int& ny, const int& nz, const int& nparticles, real_t* density){
	
	real_function_3d rho_interp;
	real_function_3d phi;
	coord_3d center;
	real_function_3d temp;
	
	set_initial_parameters(128);
	
	set_projection_precision(9, 1e-7);

	build_projected_density(world, nx, ny, nz, density, rho_interp);

	compute_potential(world, rho_interp, phi, 1e-6, 1e-8);
		
	print_potential(world, phi, 128, 128);
	
	return phi;
}

///
/// @brief Compute the gradient of the potential
///
/// @param[in]		world			Reference to the madness:world
/// @param[in]		potential		Solution of Poisson equation
/// @param[in, out]	gradient		Gradient of the solution
///

void compute_gradient(World& world, const real_function_3d& potential, vector_real_function_3d& gradient){
	
	real_derivative_3d Dx(world, 0), Dy(world, 1), Dz(world, 2);
	
	gradient[0] = Dx(potential);
	gradient[1] = Dy(potential);
	gradient[2] = Dz(potential);

}

///
/// @brief Update the position of the particles
///
/// @param[in]		world			Reference to the madness:world
/// @param[in, out]	x				x-coordinate of the particles
/// @param[in, out]	y				y-coordinate of the particles
/// @param[in, out]	z				z-coordinate of the particles
/// @param[in, out]	vx				x-component of the velocity
/// @param[in, out]	vy				y-component of the velocity
/// @param[in, out]	vz				z-component of the velocity
/// @param[in, out]	potential		Solution of Poisson equation
/// @param[in, out]	nparticles		Number of particles used
/// @param[in]		timestep		Timestep of the solver
///
/// Update position and velocity of the particles by using a Leap-Frog scheme.
///

void update_particles(World& world, real_t* x, real_t* y, real_t* z, real_t* vx, real_t* vy, real_t* vz, const int& nparticles, const real_function_3d& potential, const real_t& timestep){
		
	
	const int nx(128), ny(128), nz(128);
	
	vector_real_function_3d gradient(3);
	static const int upper_limit = nparticles;
	
	compute_gradient(world, potential, gradient);
				
	
	std::vector<real_t> grad_x, grad_y, grad_z;
	grad_x.resize(upper_limit);
	grad_y.resize(upper_limit);
	grad_z.resize(upper_limit);
	
	for(int particle = world.rank(); particle < upper_limit; particle += world.size()){

		coordT position, velocity;
		position[0] = x[particle]; position[1] = y[particle]; position[2] = z[particle];
		velocity[0] = vx[particle]; velocity[1] = vy[particle]; velocity[2] = vz[particle];
		double evaluation;


		int ndir(3);

		for(int direction(0); direction < ndir; ++direction){

			evaluation = gradient[direction].eval(position);
			velocity[direction] += evaluation * timestep;

		}

		for(int direction(0); direction < 3; ++direction){
			position[direction] += velocity[direction] * timestep;
		}

		x[particle] = position[0] + ( position[0] > nx ? -(nx-1) : (position[0] < 1 ? (nx-1) : 0) );
		y[particle] = position[1] + ( position[1] > ny ? -(ny-1) : (position[1] < 1 ? (ny-1) : 0) );
		z[particle] = position[2] + ( position[2] > nz ? -(nz-1) : (position[2] < 1 ? (nz-1) : 0) );

		vx[particle] = velocity[0]; vy[particle] = velocity[1]; vz[particle] = velocity[2];

	}

	world.gop.fence();
	
	// // Test of update using future, maybe it could help, maybe not.
	// for(int particle = world.rank(); particle < nparticles; particle += world.size()){
	//
	// 	coordT position;
	//
	// 	position[0] = x[particle]; position[1] = y[particle]; position[2] = z[particle];
	//
	// 	grad_x[particle] = gradient[0].eval(position);
	// 	grad_y[particle] = gradient[1].eval(position);
	// 	grad_z[particle] = gradient[2].eval(position);
	//
	// }
	//
	// // Force futures
	// world.gop.fence();
	// world.gop.sum(&grad_x[0], upper_limit);
	// world.gop.sum(&grad_y[0], upper_limit);
	// world.gop.sum(&grad_z[0], upper_limit);
	//
	// for(int particle = world.rank(); particle < upper_limit; particle += world.size()){
	//
	// 	coordT position, velocity;
	// 	position[0] = x[particle]; position[1] = y[particle]; position[2] = z[particle];
	// 	velocity[0] = vx[particle]; velocity[1] = vy[particle]; velocity[2] = vz[particle];
	//
	// 	velocity[0] += grad_x[particle] * timestep;
	// 	velocity[1] += grad_y[particle] * timestep;
	// 	velocity[2] += grad_z[particle] * timestep;
	//
	// 	for(int direction(0); direction < 3; ++direction){
	// 		position[direction] += velocity[direction] * timestep;
	// 	}
	//
	// 	x[particle] = position[0] + ( position[0] > nx ? -(nx-1) : (position[0] < 1 ? (nx-1) : 0) );
	// 	y[particle] = position[1] + ( position[1] > ny ? -(ny-1) : (position[1] < 1 ? (ny-1) : 0) );
	// 	z[particle] = position[2] + ( position[2] > nz ? -(nz-1) : (position[2] < 1 ? (nz-1) : 0) );
	//
	// 	vx[particle] = velocity[0]; vy[particle] = velocity[1]; vz[particle] = velocity[2];
	// }

	
}

int main(int argc, char** argv){
	
	int nx, ny, nz, nparticles, nproc;
	std::vector<real_t> x, y, z, vx, vy, vz, mass, density;
	
	if(argc == 2){
		nx = atoi(argv[1]);
		ny = nx;
		nz = nx;
	}
	
	initialize(argc, argv);
	World world(SafeMPI::COMM_WORLD);
	startup(world, argc, argv);

	if(world.rank() == 0){
		get_dim_(&nx, &ny, &nz, &nparticles, &nproc);
	}
	
	// Broadcast dimensions to all MPI ranks
	world.gop.broadcast(nx);
	world.gop.broadcast(ny);
	world.gop.broadcast(nz);
	world.gop.broadcast(nparticles);

	x.resize(nparticles);
	y.resize(nparticles);
	z.resize(nparticles);
	vx.resize(nparticles);
	vy.resize(nparticles);
	vz.resize(nparticles);
	mass.resize(nparticles);
	density.resize(nx*ny*nz);
	
	// Sanity
	if (world.rank() == 0) printf("Dimensions: %i %i %i\n", nx, ny, nz);
	if (world.rank() == 0) printf("Number of particles: %i\n", nparticles);
	if (world.rank() == 0) printf("num_procs to write file: %i\n", nproc);
	
	auto start_time = std::chrono::high_resolution_clock::now();
	
	int step(0);
	
	// Initialize density on uniform mesh 
	if(world.rank() == 0){
		part_init_(&nx, &ny, &nz, &nparticles, &nproc, &x[0], &y[0], &z[0], &vx[0], &vy[0], &vz[0], &mass[0]);
		project_density_(&nx, &ny, &nz, &nparticles, &x[0], &y[0], &z[0], &mass[0], &density[0], &step);
	}
	
	// Broadcast density to all ranks
	world.gop.broadcast(&density[0], nx*ny*nz, 0);

	world.gop.fence();
	
	real_function_3d potential = solve_potential(world, &x[0], &y[0], &z[0], nx, ny, nz, nparticles, &density[0]);
		
	// // Implementation of the complete solver, not just Poisson equation
	// nstep = 3;
	//
	// for(int step(0); step < nstep; ++step){
	//
	//
	// int step(0);
	//
	// real_function_3d potential = solve_potential(world, &x[0], &y[0], &z[0], nx, ny, nz, nparticles, &density[0]);
	//
	// update_particles(world, &x[0], &y[0], &z[0], &vx[0], &vy[0], &vz[0], nparticles, potential, timestep);
	//
	//
	// 	world.gop.fence();
	//
	// 	memset(&density[0], 0, sizeof(real_t)*nx*ny*nz);
	// }
	//
	
	world.gop.fence();
	auto overall_time = std::chrono::high_resolution_clock::now();
	if (world.rank() == 0) printf("\nOverall time: %f s\n\n", 1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(overall_time  - start_time).count());

	finalize();
	
	return 0;
}

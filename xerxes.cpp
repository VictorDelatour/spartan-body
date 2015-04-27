#include <vector>
#include <cstdio>

#include <madness/mra/mra.h>

#include "density_projector3d.hpp"

using namespace madness;

typedef double real_t;
typedef Vector<double,3> coordT;

extern "C" void get_dim_(int* nx, int* ny, int* nz, int* nparticles);
extern "C" void part_init_(const int* nx, const int* ny, const int* nz, const int* nparticles, real_t* x, real_t* y, real_t* z, real_t* vx, real_t* vy, real_t* vz, real_t* mass);
extern "C" void project_density_(const int* nx, const int* ny, const int* nz, const int* nparticles, const real_t* x, const real_t* y, const real_t* z, real_t* mass, real_t* density, const int* step);


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

void build_projected_density(World& world, const int& nx, const int& ny, const int& nz, real_t* density, real_function_3d& projected_density){
	
	real_functor_3d density_functor;
	
	density_functor = real_functor_3d(new DensityProjector(nx, ny, nz, &density[0]));
	projected_density = real_factory_3d(world).functor(density_functor);
	
}

void compute_potential(World& world, const real_function_3d& projected_density, real_function_3d& potential, const double& precision, const double& threshold){
	
	double integral, volume, mean;

	if (world.rank() == 0) printf("\tProjecting potential\n");
	
	real_convolution_3d coulomb_operator = CoulombOperator(world, precision, threshold);
	
	potential = coulomb_operator(projected_density);

	
	if (world.rank() == 0) printf("\tProjected potential\n");

	integral = (potential).trace();
	volume = FunctionDefaults<3>::get_cell_volume();
	mean = integral/volume;

	potential = potential - mean;
	if (world.rank() == 0) printf("\tNormalized\n");
	
	// if (world.rank() == 0) printf("\t#YOLO FTW");
	// real_derivative_3d Dx(world, 0), Dy(world, 1), Dz(world, 2);
	// real_function_3d deriv_x = Dx(*potential);
	
}

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

real_function_3d solve_potential(World& world, const int& nx, const int& ny, const int& nz, real_t* density){
	
	real_function_3d rho_interp;
	real_function_3d phi;
	
	if (world.rank() == 0) printf("Setup initial parameters\n");
	set_initial_parameters(nx);
	if (world.rank() == 0) printf("Set...\n\n");
	
	if (world.rank() == 0) printf("Setup projection precision\n");
	set_projection_precision(9, 1e-7);
	if (world.rank() == 0) printf("Set...\n\n");
	
	if (world.rank() == 0) printf("Build projected density\n");
	build_projected_density(world, nx, ny, nz, density, rho_interp);
	if (world.rank() == 0) printf("Built...\n\n");
	
	// if (world.rank() == 0) printf("Printing density\n");
	// print_density(world, &rho_interp, 128, nx);
	// if (world.rank() == 0) printf("Printed...\n\n");

	if (world.rank() == 0) printf("Computing potential\n");
	
	compute_potential(world, rho_interp, phi, 1e-6, 1e-8);
	
	// //
	// // potential = &phi;
	//
	// compute_potential(world, rho_interp, potential, 1e-6, 1e-8);
	if (world.rank() == 0) printf("Computed...\n\n");
	
	//
	// double temp;
	// temp = phi(5.0, 5.0, 5.0);
	// if (world.rank() == 0) printf("Eval potential %f\n", temp);
	
	// if (world.rank() == 0) printf("Printing potential\n");
	// print_potential(world, potential, 128, nx);
	// if (world.rank() == 0) printf("Printed...\n\n");	
	
	return phi;
}

void compute_gradient(World& world, const real_function_3d& potential, vector_real_function_3d& gradient){
	
	real_derivative_3d Dx(world, 0), Dy(world, 1), Dz(world, 2);
	
	// real_function_3d deriv_x = Dx(potential);
	
	// real_function_3d deriv_x = Dy(*potential);
	//
	gradient[0] = Dx(potential);
	gradient[1] = Dy(potential);
	gradient[2] = Dz(potential);
	//
}

void update_velocity(const coordT& position, coordT& velocity, const double& time_step, vector_real_function_3d& gradient){
	
	for(int direction(0); direction < 3; ++direction){
		velocity[direction] += gradient[direction].eval(position) * time_step;
	}
	
}

void update_position(coordT& position, const coordT& velocity, const double& time_step){
	
	for(int axis(0); axis < 3; ++axis){
		position[axis] += velocity[axis] * time_step;
	}
	
}

void update_particles(World& world, real_t* x, real_t* y, real_t* z, real_t* vx, real_t* vy, real_t* vz, const int& nparticles, const real_function_3d& potential, const real_t& timestep){
	
	const int nx(128), ny(128), nz(128);
	
	vector_real_function_3d gradient(3);
	// int upper_limit = 1e5;
	int upper_limit = nparticles;
	double start_time, update_time;

	if (world.rank() == 0) printf("Updating %i of %i particles...\n", upper_limit, nparticles);

	if (world.rank() == 0) printf("\tComputing gradient...\n");
	compute_gradient(world, potential, gradient);
	
	coordT temp;
	temp[0] = gradient[0](x[0], y[0], z[0]);
	temp[1] = gradient[1](x[0], y[0], z[0]);
	temp[2] = gradient[2](x[0], y[0], z[0]);
	
	

	if (world.rank() == 0) printf("\tDone.\n");

	if (world.rank() == 0) printf("\tLooping over all particles... with %i processors\n", world.size());

	start_time = wall_time();

	for(int particle = world.rank(); particle < upper_limit; particle += world.size()){

		// Is it really useful to create a coordT like this, for nothing? I don't think so...
		coordT position, velocity;
		position[0] = x[particle]; position[1] = y[particle]; position[2] = z[particle];
		velocity[0] = vx[particle]; velocity[1] = vy[particle]; velocity[2] = vz[particle];

		update_velocity(position, velocity, timestep, gradient);
		
		update_position(position, velocity, timestep);
		
		// Switch to local variable
		
		x[particle] = position[0] + ( position[0] > nx ? -(nx-1) : (position[0] < 1 ? (nx-1) : 0) ); 
		y[particle] = position[1] + ( position[1] > ny ? -(ny-1) : (position[1] < 1 ? (ny-1) : 0) ); 
		z[particle] = position[2] + ( position[2] > nz ? -(nz-1) : (position[2] < 1 ? (nz-1) : 0) );
		
		vx[particle] = velocity[0]; vy[particle] = velocity[1]; vz[particle] = velocity[2];

	}
	
	world.gop.fence();
	
	update_time = wall_time() - start_time;
	
	if (world.rank() == 0) printf("\tDone.\n\n");
	
	if (world.rank() == 0) printf("\tUpdate time was %f.\n", update_time);
		
}

int main(int argc, char** argv){
	
	int nx, ny, nz, nparticles;
	std::vector<real_t> x, y, z, vx, vy, vz, mass, density;
	double timestep;
	
	timestep = .01;
	
	initialize(argc, argv);
	World world(SafeMPI::COMM_WORLD);
	startup(world, argc, argv);

	get_dim_(&nx, &ny, &nz, &nparticles);

	x.resize(nparticles);
	y.resize(nparticles);
	z.resize(nparticles);
	vx.resize(nparticles);
	vy.resize(nparticles);
	vz.resize(nparticles);
	mass.resize(nparticles);
	density.resize(nx*ny*nz);
	

	if (world.rank() == 0) printf("Initializing particles...\n");
	part_init_(&nx, &ny, &nz, &nparticles, &x[0], &y[0], &z[0], &vx[0], &vy[0], &vz[0], &mass[0]);
	if (world.rank() == 0) printf("Done.\n");
	
	if (world.rank() == 0) printf("Particle 0 at position (%f, %f, %f) .\n", x[0], y[0], z[0]);
	
	for(int step(0); step < 5; ++step){
		
		if (world.rank() == 0) printf("Density from particles at step %i...\n", step);
		project_density_(&nx, &ny, &nz, &nparticles, &x[0], &y[0], &z[0], &mass[0], &density[0], &step);
		if (world.rank() == 0) printf("Done.\n");
	
		// potential =
		real_function_3d potential = solve_potential(world, nx, ny, nz, &density[0]);
	
	
		update_particles(world, &x[0], &y[0], &z[0], &vx[0], &vy[0], &vz[0], nparticles, potential, timestep);	
	}

	
	
	// if (world.rank() == 0) printf("Particle 0 at new position (%f, %f, %f) .\n", x[0], y[0], z[0]);
	


	finalize();
	
	return 0;
}
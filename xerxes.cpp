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

extern "C" void get_dim_(int* nx, int* ny, int* nz, int* nparticles);
extern "C" void part_init_(const int* nx, const int* ny, const int* nz, const int* nparticles, real_t* x, real_t* y, real_t* z, real_t* vx, real_t* vy, real_t* vz, real_t* mass);
extern "C" void project_density_(const int* nx, const int* ny, const int* nz, const int* nparticles, const real_t* x, const real_t* y, const real_t* z, real_t* mass, real_t* density, const int* step);

class Gaussian : public FunctionFunctorInterface<double,3> {
public:
    const coord_3d center;
    const double exponent;
    const double coefficient;
    std::vector<coord_3d> specialpt;

    Gaussian(const coord_3d& center, double exponent, double coefficient)
        : center(center), exponent(exponent), coefficient(coefficient), specialpt(1)
    {
        specialpt[0][0] = center[0];
        specialpt[0][1] = center[1];
        specialpt[0][2] = center[2];
		
		counter = new AtomicCounter();
    }
	

    // MADNESS will call this interface
    double operator()(const coord_3d& x) const {
		counter->increment();
		
		double sum = 0.0;
		for (int i(0); i<3; i++) {
			double xx = center[i]-x[i];
			xx *= 10./128.;
			sum += xx*xx;
		};

		return coefficient*exp(-exponent*sum);

    }

    // By default, adaptive projection into the spectral element basis
    // starts uniformly distributed at the initial level.  However, if
    // a function is "spiky" it may be necessary to project at a finer
    // level but doing this uniformly is expensive.  This method
    // enables us to tell MADNESS about points/areas needing deep
    // refinement (the default is no special points).
    std::vector<coord_3d> special_points() const {
        return specialpt;
    }
	
	const int get_counter() const{
		return counter->get();
	}

	
private:
	AtomicCounter* counter;
};

real_functor_3d new_gaussian(const coord_3d& origin) {

	double sigma = 1.0;

    const double exponent = 1/pow(sigma,2);
	const double coefficient = 1.0;
	
    return real_functor_3d(new Gaussian(origin, exponent, coefficient));
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

void print_density(World& world, const real_function_3d& projected_density, const int& numpts, const int& nx){
	
	const char filename_density[] = "data/spartan_density.vts";
	
	coord_3d origin;
	origin[0] = .5*(double(nx)-1) + 1;
	origin[1] = .5*(double(nx)-1) + 1;
	origin[2] = .5*(double(nx)-1) + 1;
	
	
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
	
	// if (world.rank() == 0) printf("Printed...\n\n");
	
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

void build_projected_density(World& world, const int& nx, const int& ny, const int& nz, real_t* density, real_function_3d& projected_density){
	
	if (world.rank() == 0) printf("\nDensity step\n");
	
	double access_value;
	real_functor_3d density_functor;
	
	// coord_3d center;
	// center[0] = .5 * (nx - 1.0); center[1] = .5 * (ny - 1.0); center[2] = .5 * (nz - 1.0);

	// DensityProjector numerical_gaussian(nx, ny, nz, &density[0]);
	// numerical_gaussian.reset_counter();
	
	auto start_time = std::chrono::high_resolution_clock::now();
	
	density_functor = real_functor_3d(new DensityProjector(nx, ny, nz, &density[0]));
	// density_functor = real_functor_3d( &numerical_gaussian );
	// density_functor = real_functor_3d( new Gaussian(center, 1.0, 1.0) );
	auto step_coefficients  = std::chrono::high_resolution_clock::now();
	if (world.rank() == 0) printf("\tCoefficients: %f s\n", 1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(step_coefficients - start_time).count());
	
	projected_density = real_factory_3d(world).functor(density_functor);
	auto step_projection  = std::chrono::high_resolution_clock::now();
	if (world.rank() == 0) printf("\tProjection:   %f s\n", 1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(step_projection - step_coefficients).count());

	// printf("Number of calls to numerical_gaussian: %i\n", numerical_gaussian.get_counter());

	auto printing_step  = std::chrono::high_resolution_clock::now();
	
	print_density(world, projected_density, 128, nx);
	
	auto printing_end = std::chrono::high_resolution_clock::now();
	if (world.rank() == 0) printf("\tPrinting:     %f s\n", 1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(printing_end - printing_step).count());
	
	if (world.rank() == 0) printf("Density step:         %f s\n", 1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(printing_step - start_time).count());

}

void compute_potential(World& world, const real_function_3d& projected_density, real_function_3d& potential, const double& precision, const double& threshold){
	
	if (world.rank() == 0) printf("\nPotential step\n");
	auto start_time = std::chrono::high_resolution_clock::now();
	
	double integral, volume, mean;
	
	real_convolution_3d coulomb_operator = CoulombOperator(world, precision, threshold);
	
	potential = coulomb_operator(projected_density);
	
	auto coulomb_time = std::chrono::high_resolution_clock::now();
	if (world.rank() == 0) printf("\tCoulomb step:  %f s\n", 1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(coulomb_time - start_time).count());

	// integral = potential.trace();
	// volume = FunctionDefaults<3>::get_cell_volume();
	// mean = integral/volume;
	// potential = potential - mean;
	
	auto update_time = std::chrono::high_resolution_clock::now();
	if (world.rank() == 0) printf("\tNormalization: %f s\n", 1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(update_time - coulomb_time).count());
	
	auto end_time = std::chrono::high_resolution_clock::now();
	if (world.rank() == 0) printf("Potential step:        %f s\n", 1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count());


}

real_function_3d solve_potential(World& world, real_t* x, real_t* y, real_t* z, const int& nx, const int& ny, const int& nz, const int& nparticles, real_t* density){
	
	real_function_3d rho_interp;
	real_function_3d phi;
	coord_3d center;
	real_function_3d temp;
	int limit;
	
	// if (world.rank() == 0) printf("Setup initial parameters\n");
	set_initial_parameters(nx);
	// if (world.rank() == 0) printf("Set...\n\n");
	
	// if (world.rank() == 0) printf("Setup projection precision\n");
	set_projection_precision(9, 1e-7);
	// if (world.rank() == 0) printf("Set...\n\n");
	
	// if (world.rank() == 0) printf("Build projected density\n");
	build_projected_density(world, nx, ny, nz, density, rho_interp);
	// if (world.rank() == 0) printf("Built...\n\n");
	
	// if (world.rank() == 0) printf("\tPrinting density\n");
	// print_density(world, rho_interp, 128, nx);
	// if (world.rank() == 0) printf("\tPrinted...\n\n");
	
	
	
	// if (world.rank() == 0) printf("Computing potential\n");
	
	compute_potential(world, rho_interp, phi, 1e-6, 1e-8);
		
	// if (world.rank() == 0) printf("Printing potential\n");
	print_potential(world, phi, 128, nx);
	// if (world.rank() == 0) printf("Printed...\n\n");
	
	return phi;
}

void compute_gradient(World& world, const real_function_3d& potential, vector_real_function_3d& gradient){
	
	real_derivative_3d Dx(world, 0), Dy(world, 1), Dz(world, 2);
	
	gradient[0] = Dx(potential);
	gradient[1] = Dy(potential);
	gradient[2] = Dz(potential);

}

void update_particles(World& world, real_t* x, real_t* y, real_t* z, real_t* vx, real_t* vy, real_t* vz, const int& nparticles, const real_function_3d& potential, const real_t& timestep){
		
	if (world.rank() == 0) printf("\nUpdate step\n");
	
	auto start_time = std::chrono::high_resolution_clock::now();
	
	
	const int nx(128), ny(128), nz(128);
	
	vector_real_function_3d gradient(3);
	const int upper_limit = nparticles;
	
	compute_gradient(world, potential, gradient);
	
	auto gradient_time = std::chrono::high_resolution_clock::now();
	if (world.rank() == 0) printf("\tGradient:  %f s\n",  1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(gradient_time  - start_time).count());
	
	float velocity_time(0.), position_time(0.), update_time(0.);
	//
	// Test of implementation of eval_vector
	// std::vector<real_t> x_t(1), y_t(1), z_t(1);
	// x_t[0] = x[0]; y_t[0] = y[0]; z_t[0] = z[0];

	int num_eval(1e5);

	std::vector<real_t> x_pos(num_eval), y_pos(num_eval), z_pos(num_eval);
	//
	for(int i(0); i < num_eval; ++i){
		x_pos[i] = x[i];
		y_pos[i] = y[i];
		z_pos[i] = z[i];
	}
	
	// //
	// //
	coordT pos_t;
	pos_t[0] = x_pos[0]; pos_t[1] = y_pos[0]; pos_t[2] = z_pos[0];
	// // coordT pos_t2;
	// // pos_t2[0] = x_pos[1]; pos_t2[1] = y_pos[1]; pos_t2[2] = z_pos[1];
	//
	// // Tensor<double> evaluate(num_eval);
	//
	//
	// // evaluate = gradient[0].eval_vector(x_t, y_t, z_t);
	
	
	
	Tensor<double> cell(3,2);
	std::vector<long> npoints(3);
	const int numpts(128);
	
	// if(world.rank() == 0) printf("#YOLO\n");

	for(int i(0); i < 3; ++i){
		cell(i, 0) = 1.0;
		cell(i, 1) = double(nx);
 		npoints[i] = numpts;
	}
	
	
	// if(world.rank() == 0) printf("#YOLO\n");
	
	auto fun_time = std::chrono::high_resolution_clock::now();

	Tensor<double> evaluate = gradient[0].eval_vector(x_pos, y_pos, z_pos);
	
	// Tensor<double> evaluate = gradient[0].eval_cube(cell, npoints);
	
	// evaluate = gradient[0].eval_vector(x_t, y_t, z_t);
	
	// if(world.rank() == 0) printf("#YOLO\n");
	
	auto eval_time = std::chrono::high_resolution_clock::now();
	
	double custom_time = (float)std::chrono::duration_cast<std::chrono::microseconds>(eval_time - fun_time).count();
	// if (world.rank() == 0) printf("\n Custom evaluation of %i particles: %f s\n", num_eval,  1e-6*custom_time);
	if (world.rank() == 0) printf("\n Evaluation of grad_x : %f s\n",  1e-6*custom_time);
	
	// if (world.rank() == 0) printf("\n Pos2 (%f, %f, %f)\n\n", pos_t2[0], pos_t2[1], pos_t2[2]);



	double temp = gradient[0].eval(pos_t);
	// double temp2 = gradient[0].eval(pos_t2);
	
	// if (world.rank() == 0) printf("\n Custom-built evaluation of gradient at position (%f, %f, %f): %f", x_pos[0], 	y_pos[0], 	z_pos[0], 	evaluate[0]);
	// if (world.rank() == 0) printf("\n Standard evaluation of gradient at position (%f, %f, %f): %f\n\n", pos_t[0], 	pos_t[1], 	pos_t[2], 	temp);
	
	
	// if (world.rank() == 0) printf("\n Custom-built evaluation of gradient at position (%f, %f, %f): %f", x_pos[1], 	y_pos[1], 	z_pos[1], 	evaluate[1]);
	// if (world.rank() == 0) printf("\n Standard evaluation of gradient at position (%f, %f, %f): %f\n\n", pos_t2[0], pos_t2[1], 	pos_t2[2], 	temp2);
		
	// for(int particle = world.rank(); particle < upper_limit; particle += world.size()){
	//
	// 	coordT position, velocity;
	// 	position[0] = x[particle]; position[1] = y[particle]; position[2] = z[particle];
	// 	velocity[0] = vx[particle]; velocity[1] = vy[particle]; velocity[2] = vz[particle];
	//
	// 	auto vel_time = std::chrono::high_resolution_clock::now();
	// 	int ndir(3);
	// 	// Potentially, evaluating 3 functions could be problematic
	// 	// Based on the loading/unloading of coefficients
	// 	for(int direction(0); direction < ndir; ++direction){
	// 		velocity[direction] += gradient[direction].eval(position) * timestep;
	// 	}
	//
	// 	auto pos_time = std::chrono::high_resolution_clock::now();
	// 	for(int direction(0); direction < 3; ++direction){
	// 		position[direction] += velocity[direction] * timestep;
	// 	}
	//
	// 	auto up_time = std::chrono::high_resolution_clock::now();
	//
	// 	x[particle] = position[0] + ( position[0] > nx ? -(nx-1) : (position[0] < 1 ? (nx-1) : 0) );
	// 	y[particle] = position[1] + ( position[1] > ny ? -(ny-1) : (position[1] < 1 ? (ny-1) : 0) );
	// 	z[particle] = position[2] + ( position[2] > nz ? -(nz-1) : (position[2] < 1 ? (nz-1) : 0) );
	//
	// 	vx[particle] = velocity[0]; vy[particle] = velocity[1]; vz[particle] = velocity[2];
	//
	// 	auto final_time = std::chrono::high_resolution_clock::now();
	//
	// 	velocity_time += (float)std::chrono::duration_cast<std::chrono::microseconds>(pos_time  - vel_time).count();
	// 	position_time += (float)std::chrono::duration_cast<std::chrono::microseconds>(up_time - pos_time).count();
	// 	update_time += (float)std::chrono::duration_cast<std::chrono::microseconds>(final_time - up_time).count();
	//
	// }
	
	if (world.rank() == 0) printf("\tVelocity:  %f s\n",  1e-6*velocity_time );
	if (world.rank() == 0) printf("\tPosition:  %f s\n",  1e-6*position_time );
	if (world.rank() == 0) printf("\tUpdate:    %f s\n",  1e-6*update_time );
	
	auto end_time = std::chrono::high_resolution_clock::now();
	if (world.rank() == 0) printf("Update step: %f s\n",  1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(end_time  - start_time).count());
			
}

int main(int argc, char** argv){
	
	int nx, ny, nz, nparticles;
	int nstep;
	std::vector<real_t> x, y, z, vx, vy, vz, mass, density;
	double timestep;
	
	// Put algorithm to do adaptive timestepping
	timestep = 5;
	
	if(argc == 3){
		nparticles = atoi(argv[1]);
		nx = atoi(argv[2]);
		ny = nx;
		nz = nx;
	}
	
	initialize(argc, argv);
	World world(SafeMPI::COMM_WORLD);
	startup(world, argc, argv);

	// get_dim_(&nx, &ny, &nz, &nparticles);

	x.resize(nparticles);
	y.resize(nparticles);
	z.resize(nparticles);
	vx.resize(nparticles);
	vy.resize(nparticles);
	vz.resize(nparticles);
	mass.resize(nparticles);
	density.resize(nx*ny*nz);
	
	if (world.rank() == 0) printf("Dimensions: %i %i %i\n", nx, ny, nz);
	if (world.rank() == 0) printf("Number of particles: %i\n", nparticles);
	
	// start_time = wall_time();
	auto start_time = std::chrono::high_resolution_clock::now();
	
	part_init_(&nx, &ny, &nz, &nparticles, &x[0], &y[0], &z[0], &vx[0], &vy[0], &vz[0], &mass[0]);
	
	
	world.gop.fence();
	// auto init_time = std::chrono::high_resolution_clock::now();
	// if (world.rank() == 0) printf("\nInitialization time: %f s\n\n", 1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(init_time - start_time).count());
	
	nstep = 3;
	
	// for(int step(0); step < nstep; ++step){
	//
	// 	world.gop.fence();
		// auto step_start_time = std::chrono::high_resolution_clock::now();
	//
		int step(0);
	// 	//
		project_density_(&nx, &ny, &nz, &nparticles, &x[0], &y[0], &z[0], &mass[0], &density[0], &step);
	// 	//
	//
	// 	world.gop.fence();
		// auto step_density_time = std::chrono::high_resolution_clock::now();
		// if (world.rank() == 0) printf("\tDensity %i: %f s\n", step, 1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(step_density_time - step_start_time).count());
	//
	// 	//
	// 	// real_function_3d potential = solve_potential(world, nx, ny, nz, &density[0]);
		real_function_3d potential = solve_potential(world, &x[0], &y[0], &z[0], nx, ny, nz, nparticles, &density[0]);
	// 	//
	//
		// world.gop.fence();
		// auto potential_time  = std::chrono::high_resolution_clock::now();
		// if (world.rank() == 0) printf("\tPotential %i: %f s\n", step, 1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(step_potential_time  - step_density_time).count());
		// if (world.rank() == 0) printf("\nFinished\n\n");
	//
	// 	//
		update_particles(world, &x[0], &y[0], &z[0], &vx[0], &vy[0], &vz[0], nparticles, potential, timestep);
	// 	//
	//
	// 	world.gop.fence();
		// auto update_time = std::chrono::high_resolution_clock::now();
		// if (world.rank() == 0) printf("\tUpdate %i: %f s\n", step,  1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(update_time  - potential_time).count());
	//
	// 	//
	// 	memset(&density[0], 0, sizeof(real_t)*nx*ny*nz);
	// 	//
	//
	// 	world.gop.fence();
	// 	auto step_finish_time = std::chrono::high_resolution_clock::now();
	// 	if (world.rank() == 0) printf("\nStep %i took %f s\n\n", step, 1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(step_finish_time  - step_start_time).count());
	//
	// }
	//
	world.gop.fence();
	auto overall_time = std::chrono::high_resolution_clock::now();
	if (world.rank() == 0) printf("\nOverall time: %f s\n\n", 1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(overall_time  - start_time).count());

	finalize();
	
	printf("#YOLO");
	
	return 0;
}

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

extern "C" void get_dim_(int* nx, int* ny, int* nz, int* nparticles, int* nproc);
extern "C" void part_init_(const int* nx, const int* ny, const int* nz, const int* nparticles, const int* nproc, real_t* x, real_t* y, real_t* z, real_t* vx, real_t* vy, real_t* vz, real_t* mass);
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
		// counter->increment();
		
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

const double p( 1/( 4 * constants::pi ) );//*6.67384e-11));

double gaussian_density(const coord_3d& r){
	double sum = .0;
	for (int i(0); i < 3; ++i){
		double xx = r[i];
		sum += xx*xx;
	}
	
	return -(4*sum - 6) * p * exp(-sum);
}

real_functor_3d new_gaussian(const coord_3d& origin) {

	double sigma = 1.0;

    const double exponent = 1/pow(sigma,2);
	const double coefficient = 1.0;
	
    return real_functor_3d(new Gaussian(origin, exponent, coefficient));
}

void set_initial_parameters(const int& nx){
	
	BoundaryConditions<3> bc(BC_PERIODIC);
	
	FunctionDefaults<3>::set_cubic_cell(1.0, static_cast<double>(nx));
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

void initialize_mod_density(const int& nx, const int& ny, const int& nz, real_t* density){
	
	coordT position;
	
	// double left(-5.), right(5.);
	double left(1.0), right(static_cast<double>(nx));
	double length(right-left);
	

	// position[2] = -.5*(nz + 1.);
	position[2] = left;
	
	real_t px(length / static_cast<real_t>(nx));
	real_t py(length / static_cast<real_t>(ny));
	real_t pz(length / static_cast<real_t>(nz));
	
	int index_y, index_z;
	
	for(int idz(0); idz < nz; ++idz){

		index_z = idz * nx * ny;

		// position[1] = -.5*(ny + 1.);
		position[1] = left;
		position[2] += px;

		for(int idy(0); idy < ny; ++ idy){

			index_y = index_z + idy * nx;

			// position[0] = -.5*(nx + 1.);
			position[0] = left;
			position[1] += py;

			for(int idx(0); idx < nx; ++idx){

				// Position ranges from 1 to 128
				position[0] += pz;

				density[index_y + idx] = 1e5*gaussian_density(position);
				// printf("Density (%f, %f, %f): %f\n", position[0], position[1], position[2], density[index_y + idx]);
			}
		}
	}
}

void build_projected_density(World& world, const int& nx, const int& ny, const int& nz, real_t* density, real_function_3d& projected_density){
	
	if (world.rank() == 0) printf("\nDensity step\n");

	auto start_time = std::chrono::high_resolution_clock::now();
		
	real_functor_3d density_functor = real_functor_3d(new DensityProjector(nx, ny, nz, &density[0]));
	
	auto step_coefficients  = std::chrono::high_resolution_clock::now();
	if (world.rank() == 0) printf("\tCoefficients: %f s\n", 1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(step_coefficients - start_time).count());
	
	projected_density = real_factory_3d(world).functor(density_functor);
	auto step_projection  = std::chrono::high_resolution_clock::now();
	if (world.rank() == 0) printf("\tProjection:   %f s\n", 1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(step_projection - step_coefficients).count());

	auto printing_step  = std::chrono::high_resolution_clock::now();
	
	// print_density(world, projected_density, 128, nx);
	print_density(world, projected_density, 128, 128);
	
	auto printing_end = std::chrono::high_resolution_clock::now();
	if (world.rank() == 0) printf("\tPrinting:     %f s\n", 1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(printing_end - printing_step).count());
	
	if (world.rank() == 0) printf("Density step:         %f s\n", 1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(printing_step - start_time).count());

}

void compute_potential(World& world, const real_function_3d& projected_density, real_function_3d& potential, const double& precision, const double& threshold){
	
	if (world.rank() == 0) printf("\nPotential step\n");
	auto start_time = std::chrono::high_resolution_clock::now();
	
	// double integral, volume, mean;
	
	real_convolution_3d coulomb_operator = CoulombOperator(world, precision, threshold);
	
	potential = coulomb_operator(projected_density);
	
	auto coulomb_time = std::chrono::high_resolution_clock::now();
	if (world.rank() == 0) printf("\tCoulomb step:  %f s\n", 1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(coulomb_time - start_time).count());

	// integral = potential.trace();
	// integral = potential.norm2();
	// volume = FunctionDefaults<3>::get_cell_volume();
	// mean = integral*integral/volume;
	//
	// if (world.rank() == 0) printf("Integral = %f\n", integral);
	// if (world.rank() == 0) printf("Volume   = %f\n", volume);
	// if (world.rank() == 0) printf("Mean		= %f\n", mean);
	//
	
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
	
	set_initial_parameters(128);
	
	set_projection_precision(9, 1e-7);

	build_projected_density(world, nx, ny, nz, density, rho_interp);

	compute_potential(world, rho_interp, phi, 1e-6, 1e-8);
		
	// print_potential(world, phi, 128, nx);
	print_potential(world, phi, 128, 128);
	
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
	static const int upper_limit = nparticles;
	
	compute_gradient(world, potential, gradient);
	
	auto gradient_time = std::chrono::high_resolution_clock::now();
	if (world.rank() == 0) printf("\tGradient:  %f s\n",  1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(gradient_time  - start_time).count());
	
	// float velocity_time(0.), position_time(0.), update_time(0.);
			
	
	auto eval_time = std::chrono::high_resolution_clock::now();
	
	std::vector<real_t> grad_x, grad_y, grad_z;
	grad_x.resize(upper_limit);
	grad_y.resize(upper_limit);
	grad_z.resize(upper_limit);
	
	if (world.rank() == 0) printf(" Computing futures\n");
	// Use futures to speed up stuff?
	
	for(int particle = world.rank(); particle < nparticles; particle += world.size()){
		
		coordT position;
		
		position[0] = x[particle]; position[1] = y[particle]; position[2] = z[particle];
		
		grad_x[particle] = gradient[0].eval(position);
		// grad_y[particle] = gradient[1].eval(position);
		// grad_z[particle] = gradient[2].eval(position);
		// printf("Particle %i\n", particle);
		
	}
	
	if (world.rank() == 0) printf(" Forcing futures\n");
	
	// Force futures
	world.gop.fence();
	world.gop.sum(&grad_x[0], upper_limit);
	// world.gop.sum(&grad_y[0], upper_limit);
	// world.gop.sum(&grad_z[0], upper_limit);
	
	auto up_time = std::chrono::high_resolution_clock::now();
	
	if (world.rank() == 0) printf("Evaluation time: %f s\n",  1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(up_time - eval_time).count());
	
	for(int particle = world.rank(); particle < upper_limit; particle += world.size()){
		
		coordT position, velocity;
		position[0] = x[particle]; position[1] = y[particle]; position[2] = z[particle];
		velocity[0] = vx[particle]; velocity[1] = vy[particle]; velocity[2] = vz[particle];
		
		velocity[0] += grad_x[particle] * timestep;
		velocity[1] += grad_y[particle] * timestep;
		velocity[2] += grad_z[particle] * timestep;
	
		for(int direction(0); direction < 3; ++direction){
			position[direction] += velocity[direction] * timestep;
		}

		x[particle] = position[0] + ( position[0] > nx ? -(nx-1) : (position[0] < 1 ? (nx-1) : 0) );
		y[particle] = position[1] + ( position[1] > ny ? -(ny-1) : (position[1] < 1 ? (ny-1) : 0) );
		z[particle] = position[2] + ( position[2] > nz ? -(nz-1) : (position[2] < 1 ? (nz-1) : 0) );

		vx[particle] = velocity[0]; vy[particle] = velocity[1]; vz[particle] = velocity[2];
	}
	
	// for(int particle = world.rank(); particle < upper_limit; particle += world.size()){
//
// 		coordT position, velocity;
// 		position[0] = x[particle]; position[1] = y[particle]; position[2] = z[particle];
// 		velocity[0] = vx[particle]; velocity[1] = vy[particle]; velocity[2] = vz[particle];
// 		double evaluation;
//
// 		// auto vel_time = std::chrono::high_resolution_clock::now();
// 		int ndir(3);
// 		// Potentially, evaluating 3 functions could be problematic
// 		// Based on the loading/unloading of coefficients
// 		for(int direction(0); direction < ndir; ++direction){
//
// 			evaluation = gradient[direction].eval(position);
// 			velocity[direction] += evaluation * timestep;
//
// 		}
//
// 		// auto pos_time = std::chrono::high_resolution_clock::now();
// 		for(int direction(0); direction < 3; ++direction){
// 			position[direction] += velocity[direction] * timestep;
// 		}
//
// 		// auto up_time = std::chrono::high_resolution_clock::now();
//
// 		x[particle] = position[0] + ( position[0] > nx ? -(nx-1) : (position[0] < 1 ? (nx-1) : 0) );
// 		y[particle] = position[1] + ( position[1] > ny ? -(ny-1) : (position[1] < 1 ? (ny-1) : 0) );
// 		z[particle] = position[2] + ( position[2] > nz ? -(nz-1) : (position[2] < 1 ? (nz-1) : 0) );
//
// 		vx[particle] = velocity[0]; vy[particle] = velocity[1]; vz[particle] = velocity[2];
//
// 		// auto final_time = std::chrono::high_resolution_clock::now();
//
// 		// velocity_time += (float)std::chrono::duration_cast<std::chrono::microseconds>(pos_time  - vel_time).count();
// 		// position_time += (float)std::chrono::duration_cast<std::chrono::microseconds>(up_time - pos_time).count();
// 		// update_time += (float)std::chrono::duration_cast<std::chrono::microseconds>(final_time - up_time).count();
//
// 		// printf("Particle %i\n", particle);
// 	}
//
	world.gop.fence();
	
	// if (world.rank() == 0) printf("\tVelocity:  %f s\n",  1e-6*velocity_time );
	// if (world.rank() == 0) printf("\tPosition:  %f s\n",  1e-6*position_time );
	// if (world.rank() == 0) printf("\tUpdate:    %f s\n",  1e-6*update_time );
	
	auto end_time = std::chrono::high_resolution_clock::now();
	if (world.rank() == 0) printf("Update step: %f s\n",  1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(end_time  - start_time).count());
			
}

int main(int argc, char** argv){
	
	int nx, ny, nz, nparticles, nproc;
	// int nstep;
	std::vector<real_t> x, y, z, vx, vy, vz, mass, density;
	
	if(argc == 2){
		nx = atoi(argv[1]);
		ny = nx;
		nz = nx;
	}
	
	initialize(argc, argv);
	World world(SafeMPI::COMM_WORLD);
	startup(world, argc, argv);

	// Better to only read once all the data
	if(world.rank() == 0){
		get_dim_(&nx, &ny, &nz, &nparticles, &nproc);
	}
	
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
	
	
	if (world.rank() == 0) printf("Dimensions: %i %i %i\n", nx, ny, nz);
	if (world.rank() == 0) printf("Number of particles: %i\n", nparticles);
	if (world.rank() == 0) printf("num_procs to write file: %i\n", nproc);
	
	auto start_time = std::chrono::high_resolution_clock::now();
	
	int step(0);
	
	if(world.rank() == 0){
		part_init_(&nx, &ny, &nz, &nparticles, &nproc, &x[0], &y[0], &z[0], &vx[0], &vy[0], &vz[0], &mass[0]);
		project_density_(&nx, &ny, &nz, &nparticles, &x[0], &y[0], &z[0], &mass[0], &density[0], &step);
	}
	
	world.gop.broadcast(&density[0], nx*ny*nz, 0);

	world.gop.fence();
		
	// nstep = 3;
	//
	// for(int step(0); step < nstep; ++step){
	//
	//
	// int step(0);
	//
	real_function_3d potential = solve_potential(world, &x[0], &y[0], &z[0], nx, ny, nz, nparticles, &density[0]);
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

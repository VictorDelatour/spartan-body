#include <vector>
#include <cstdio>

#include <madness/mra/mra.h>

#include "solve_potential.hpp"
#include "update_particles.hpp"

typedef double real_t;

extern "C" void get_dim_(int* nx, int* ny, int* nz, int* nparticles);
extern "C" void part_init_(const int* nx, const int* ny, const int* nz, const int* nparticles, real_t* x, real_t* y, real_t* z, real_t* vx, real_t* vy, real_t* vz, real_t* mass);
extern "C" void project_density_(const int* nx, const int* ny, const int* nz, const int* nparticles, const real_t* x, const real_t* y, const real_t* z, real_t* mass, real_t* density);

using namespace madness;

int main(int argc, char** argv){
	
	int nx, ny, nz, nparticles;
	std::vector<real_t> x, y, z, vx, vy, vz, mass, density;
	real_function_3d potential;
	double timestep;
	
	timestep = .5;
	
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

	if (world.rank() == 0) printf("Density from particles...\n");
	project_density_(&nx, &ny, &nz, &nparticles, &x[0], &y[0], &z[0], &mass[0], &density[0]);
	if (world.rank() == 0) printf("Done.\n");

	solve_potential(world, nx, ny, nz, &density[0], &potential);

	// if (world.rank() == 0) printf("Potential evaluation: %f\n", potential(5.0, 5.0, 5.0));

	// update_particles(&world, &x[0], &y[0], &z[0], &vx[0], &vy[0], &vz[0], nparticles, potential, timestep);

	if (world.rank() == 0) printf("Particles updated, WHAT ARE YOU WAITING FOR!\n");

	finalize();
	
	return 0;
}
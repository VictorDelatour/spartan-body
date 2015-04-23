#include <vector>
#include <cstdio>

#include <madness/mra/mra.h>

#include "solve_potential.hpp"

typedef double real_t;

extern "C" void get_dim_(int* nx, int* ny, int* nz, int* nparticles);
extern "C" void part_init_(const int* nx, const int* ny, const int* nz, const int* nparticles, real_t* x, real_t* y, real_t* z, real_t* mass);
extern "C" void project_density_(const int* nx, const int* ny, const int* nz, const int* nparticles, const real_t* x, const real_t* y, const real_t* z, real_t* mass, real_t* density);

using namespace madness;

int main(int argc, char** argv){
	
	int nx, ny, nz, nparticles;
	std::vector<real_t> x, y, z, mass, density;
	
	
	initialize(argc, argv);
	World world(SafeMPI::COMM_WORLD);
	startup(world, argc, argv);
		
	get_dim_(&nx, &ny, &nz, &nparticles);
	
	x.resize(nparticles);
	y.resize(nparticles);
	z.resize(nparticles);
	mass.resize(nparticles);
	density.resize(nx*ny*nz);
		
	part_init_(&nx, &ny, &nz, &nparticles, &x[0], &y[0], &z[0], &mass[0]);
	
	project_density_(&nx, &ny, &nz, &nparticles, &x[0], &y[0], &z[0], &mass[0], &density[0]);
	
	solve_potential(&world, nx, ny, nz, &density[0]);
	
	finalize();
	
	return 0;
}
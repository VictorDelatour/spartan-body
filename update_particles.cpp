#include "update_particles.hpp"

void update_particles(World* world, real_t* x, real_t* y, real_t* z, const int& nparticles, const real_function_3d& potential){
	
	real_function_3d dx_phi, dy_phi, dz_phi;
	
}

void compute_gradient(World* world, const real_function_3d& potential, real_function_3d& dx_potential, real_function_3d& dy_potential, real_function_3d& dz_potential){
	
	real_derivative_3d Dx(*world, Direction::X), Dy(*world, Direction::Y),Dz(*world, Direction::Z);
	
	dx_potential = Dx(potential);
	dy_potential = Dy(potential);
	dz_potential = Dz(potential);
	
}
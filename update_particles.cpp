#include "update_particles.hpp"

void update_particles(World* world, real_t* x, real_t* y, real_t* z, real_t* vx, real_t* vy, real_t* vz, const int& nparticles, const real_function_3d& potential, const real_t& timestep){
	
	vector_real_function_3d gradient(3);
	// real_function_3d dx_phi, dy_phi, dz_phi;
	
	compute_gradient(world, potential, gradient);
	// compute_gradient(*world, potential, dx_phi, dy_phi, dz_phi);
	
	for(int particle(0); particle < nparticles; ++particle){
		
		update_velocity(&x[particle], &y[particle], &z[particle], &vx[particle], &vy[particle], &vz[particle], timestep, gradient);
		
		update_position(&x[particle], &y[particle], &z[particle], &vx[particle], &vy[particle], &vz[particle], timestep);
		
	}
	
}

void compute_gradient(World* world, const real_function_3d& potential, vector_real_function_3d& gradient){
	
	real_derivative_3d Dx(*world, 0), Dy(*world, 1), Dz(*world, 2);
	
	gradient[0] = Dx(potential);
	gradient[1] = Dy(potential);
	gradient[2] = Dz(potential);
	
}

// void compute_gradient(World* world, const real_function_3d& potential, real_function_3d& dx_potential, real_function_3d& dy_potential, real_function_3d& dz_potential){
//
// 	real_derivative_3d Dx(*world, 0), Dy(*world, 1),Dz(*world, 2);
//
// 	dx_potential = Dx(potential);
// 	dy_potential = Dy(potential);
// 	dz_potential = Dz(potential);
//
// }

// void update_velocity(real_t* vx, real_t* vy, real_t* vz, const double& time_step, const real_function_3d& dx_potential, const real_function_3d& dy_potential, const real_function_3d& dz_potential){
//
// }

void update_velocity(real_t* x, real_t* y, real_t* z, real_t* vx, real_t* vy, real_t* vz, const real_t& time_step, vector_real_function_3d& gradient){
	
	 *vx += gradient[0](*x, *y, *z) * time_step;
	 *vy += gradient[1](*x, *y, *z) * time_step;
	 *vz += gradient[2](*x, *y, *z) * time_step;

}

void update_position(real_t* x, real_t* y, real_t* z, real_t* vx, real_t* vy, real_t* vz, const real_t& time_step){
	
	*x += (*vx) * time_step;
	*y += (*vy) * time_step;
	*z += (*vz) * time_step;
	
}
#include <madness/mra/mra.h>
// #include <vector>

using namespace madness;

typedef double real_t;
typedef Vector<double,3> coordT;

enum class Direction {X = 0, Y = 1, Z = 2};



void update_particles(World* world, real_t* x, real_t* y, real_t* z, real_t* vx, real_t* vy, real_t* vz, const int& nparticles, const real_function_3d& potential, const double& timestep);
// void compute_gradient(World* world, const real_function_3d& potential, real_function_3d& dx_potential, real_function_3d& dy_potential, real_function_3d& dz_potential);
void compute_gradient(World* world, const real_function_3d& potential, vector_real_function_3d& gradient);
// void update_velocity(real_t* vx, real_t* vy, real_t* vz, const double& time_step, const real_function_3d& dx_potential, const real_function_3d& dy_potential, const real_function_3d& dz_potential);

// void update_velocity(real_t* x, real_t* y, real_t* z, real_t* vx, real_t* vy, real_t* vz, const double& time_step, vector_real_function_3d& gradient);
// void update_position(real_t* x, real_t* y, real_t* z, real_t* vx, real_t* vy, real_t* vz, const double& time_step);

void update_velocity(const coordT& position, coordT& velocity, const double& time_step, vector_real_function_3d& gradient);
void update_position(coordT& position, const coordT& velocity, const double& time_step);
#include <madness/mra/mra.h>

typedef double real_t;

enum class Direction {X = 0, Y = 1, Z = 2};

using namespace madness;

void update_particles(World* world, real_t* x, real_t* y, real_t* z, const int& nparticles, const real_function_3d& potential);
void compute_gradient(World* world, const real_function_3d& potential, real_function_3d& dx_potential, real_function_3d& dy_potential, real_function_3d& dz_potential);
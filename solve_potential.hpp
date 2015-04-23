#include <vector>
#include "density_projector3d.hpp"
#include <madness/mra/mra.h>

typedef double real_t;

using namespace madness;

void solve_potential(World *world, const int& nx, const int& ny, const int& nz, real_t* density);
void set_initial_parameters(const int& nx);
void set_projection_precision(const int& order, const double& threshold);
void build_projected_density(World *world, const int& nx, const int& ny, const int& nz, real_t* density, real_function_3d& projected_density);
void print_density(World* world, const real_function_3d& projected_density, const int& numpts, const int& nx);
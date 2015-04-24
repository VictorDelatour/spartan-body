#include <vector>

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

#include "density_projector3d.hpp"

typedef double real_t;

using namespace madness;

void solve_potential(World& world, const int& nx, const int& ny, const int& nz, real_t* density, real_function_3d& potential);
void set_initial_parameters(const int& nx);
void set_projection_precision(const int& order, const double& threshold);
void build_projected_density(World& world, const int& nx, const int& ny, const int& nz, real_t* density, real_function_3d& projected_density);
void compute_potential(World& world, real_function_3d& projected_density, real_function_3d& potential, const double& precision, const double& threshold);
void print_density(World& world, const real_function_3d& projected_density, const int& numpts, const int& nx);
void print_potential(World& world, const real_function_3d& potential, const int& numpts, const int& nx);
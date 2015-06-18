#include "density_projector3d.hpp"

DensityProjector::DensityProjector(const int nx, const int ny, const int nz, double* data):
nx(nx), ny(ny), nz(nz), data(data){
	
	counter = new AtomicCounter();
	
	// Filtering along rows, the space between two consecutive elements is simply 1 (adjacent)
	inplace_filtering(Axis::X);
	// Filtering along columns, the space between two consecutive elements is nx
	inplace_filtering(Axis::Y);
	// Filtering along slabs, the space between two consecutive elements is the nx*ny plane
	inplace_filtering(Axis::Z);
	
	
}

DensityProjector::~DensityProjector(){
	// get_counter();
	// printf("Number of accesses %i\n", counter->get());
	delete data;
	delete counter;
}

const int DensityProjector::get_counter() const{
	return counter->get();
}
void DensityProjector::reset_counter(){
	counter->reset();
}

const double DensityProjector::test_performance(const madness::coord_3d& x, const int& npoints) const{
	
	double 	returnvalue;
	double  interx, intery;			// Temporary variable to store the interpolated values along x and y axis
	int 	disp;					// Position

	double xweights[4], yweights[4], zweights[4];
	int xpositions[4], ypositions[4], zpositions[4];
	
	auto start_time = std::chrono::high_resolution_clock::now();
	
	for(int iter(0); iter < npoints; ++iter){
		
		get_weights_and_position(Axis::X, x, &xweights[0], &xpositions[0]);
		get_weights_and_position(Axis::Y, x, &yweights[0], &ypositions[0]);
		get_weights_and_position(Axis::Z, x, &zweights[0], &zpositions[0]);
	}
	
	auto weight_time = std::chrono::high_resolution_clock::now();
	
	for(int iter(0); iter < npoints; ++iter){
		
		returnvalue = 0.0;

		for(int k(0); k <= 3; ++k){

			intery = 0.0; 

			for(int j(0); j <= 3; ++j){

				disp = (ypositions[j] + zpositions[k] * ny) * nx;
				interx = xweights[0] * data[disp + xpositions[0]] + xweights[1] * data[disp + xpositions[1]] + xweights[2] * data[disp + xpositions[2]] + xweights[3] * data[disp + xpositions[3]];
				
				intery += yweights[j] * interx;
			}
			returnvalue += zweights[k] * intery;
		}
	}
	
	auto end_time = std::chrono::high_resolution_clock::now();
	
	
	printf("Return value = %f \n", returnvalue);
	double weights_time = 1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(weight_time - start_time).count();
	double interp_time = 1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(end_time - weight_time).count();
	double total_time = 1e-3*(float)std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	printf("Weights, %i points: %f s\n", npoints, weights_time);
	printf("Linear interpolation, %i points: %f s\n", npoints, interp_time);
	printf("Total time, %i points: %f s\n", npoints, total_time);
	
	return total_time;
}

double DensityProjector::operator()(const madness::coord_3d& x) const{
	
	// counter->increment();

	double 	returnvalue;
	double  interx, intery;			// Temporary variable to store the interpolated values along x and y axis
	int 	disp;					// Position
	
	double xweights[4], yweights[4], zweights[4];
	int xpositions[4], ypositions[4], zpositions[4];

	get_weights_and_position(Axis::X, x, &xweights[0], &xpositions[0]);
	get_weights_and_position(Axis::Y, x, &yweights[0], &ypositions[0]);
	get_weights_and_position(Axis::Z, x, &zweights[0], &zpositions[0]);


	returnvalue = 0.0;

	// Actually compute the value
	for(int k(0); k <= 3; ++k){

		intery = 0.0; // Interpolated value along x-axis

		for(int j(0); j <= 3; ++j){

			interx = 0.0; // Interpolated value along x-axis
			//
			// for(int i(0); i <= 3; ++i){
			// 	disp = xpositions[i] + (ypositions[j] + zpositions[k] * ny) * nx; // Position
			//
			// 	interx += xweights[i] * data[disp];
			// }
			
			// Loop unrolling
			disp = (ypositions[j] + zpositions[k] * ny) * nx;
			interx = xweights[0] * data[disp + xpositions[0]] + xweights[1] * data[disp + xpositions[1]] + xweights[2] * data[disp + xpositions[2]] + xweights[3] * data[disp + xpositions[3]];

			intery += yweights[j] * interx;
		}
		returnvalue += zweights[k] * intery;
	}

	return returnvalue;
}

int DensityProjector::inplace_filtering(Axis axis){
	
	int size_dim1, size_dim2, shift_dim1, shift_dim2, process_size, access_shift;
	
	switch(axis){
		case Axis::X: {
			size_dim1 = ny;
			size_dim2 = nz;
			shift_dim1 = nx;
			shift_dim2 = nx*ny;
			process_size = nx;
			access_shift = 1;
			break;
		}
		case Axis::Y:{
			size_dim1 = nx;
			size_dim2 = nz;
			shift_dim1 = 1;
			shift_dim2 = nx*ny;
			process_size = ny;
			access_shift = nx;
			break;
		}
		case Axis::Z:{
			size_dim1 = nx;
			size_dim2 = ny;
			shift_dim1 = 1;
			shift_dim2 = nx;
			process_size = nz;
			access_shift = nx*ny;
			break;
		}
	}
	
	
	for(int i(0); i < size_dim1; ++i){
		for(int j(0); j < size_dim2; ++j){
			get_coef(&data[i * shift_dim1 + j * shift_dim2], process_size, access_shift);
		}
	}

	return 0;
	
}

int DensityProjector::gaussian_filtering(Axis axis){
	
	int size_dim1, size_dim2, shift_dim1, shift_dim2, process_size, access_shift;
	int ij_shift;
	
	switch(axis){
		case Axis::X: {
			size_dim1 = ny;
			size_dim2 = nz;
			shift_dim1 = nx;
			shift_dim2 = nx*ny;
			process_size = nx;
			access_shift = 1;
			break;
		}
		case Axis::Y:{
			size_dim1 = nx;
			size_dim2 = nz;
			shift_dim1 = 1;
			shift_dim2 = nx*ny;
			process_size = ny;
			access_shift = nx;
			break;
		}
		case Axis::Z:{
			size_dim1 = nx;
			size_dim2 = ny;
			shift_dim1 = 1;
			shift_dim2 = nx;
			process_size = nz;
			access_shift = nx*ny;
			break;
		}
	}
	
	real_t first, last, prev, temp, p;
	p = .25;
	
	for(int i(0); i < size_dim1; ++i){
		for(int j(0); j < size_dim2; ++j){
			ij_shift = i * shift_dim1 + j * shift_dim2;
			
			first = data[ij_shift];
			last = data[ij_shift + access_shift * (process_size-1)];
			
			prev = p * (last + 2*first + data[ij_shift + access_shift]);
			
			for(int k(1); k < process_size - 1; ++k ){
				temp = p * (data[ij_shift + (k-1)*access_shift] + 2*data[ij_shift + k*access_shift] + data[ij_shift + (k+1)*access_shift]);
				data[ij_shift + (k-1)*access_shift] = prev;
				prev = temp;
			}
			
			data[ij_shift + access_shift * (process_size-1)] = p * (data[ij_shift + access_shift * (process_size-2)] + 2*data[ij_shift + access_shift * (process_size-1)] + first);
			
		}
	}
	
	return 0;
	
}

int DensityProjector::get_weights_and_position(Axis axis, const madness::coord_3d& x, double* weights, int* position) const{
	
	int pos, index, numel, temp_position;
	double alpha;
	
	numel = ( (axis == Axis::X) ? nx : ( (axis == Axis::Y) ? ny : nz) );
	pos = ( (axis == Axis::X) ? 0 : ( (axis == Axis::Y) ? 1 : 2) );
	
	// Is this slow?
	// index = (int) floor(x[pos] - 1);
	index = static_cast<int> (x[pos] - 1.); // Faster, plus you should be safe, as x belongs to [1, 128]
	
	alpha = x[pos] - (index + 1.); 	
	
	double ialp(1. - alpha);
	double ialp2(ialp * ialp);
	double alp2(alpha * alpha);

	weights[0] = 1./6. * ialp * ialp2;
	weights[1] = 2./3. - .5 * alp2 * (2. - alpha);
	weights[2] = 2./3. - .5 * ialp2 * (1. + alpha);
	weights[3] = 1./6. * alp2 * alpha;

	// weights[0] = 1./6. * (1. - alpha) * (1. - alpha) * (1. - alpha);
	// weights[1] = 2./3. - .5 * alpha * alpha * (2. - alpha);
	// weights[2] = 2./3. - .5 * (1. - alpha) * (1. - alpha) * (1. + alpha);
	// weights[3] = 1./6. * alpha * alpha * alpha;
	
	for(int k(0); k <= 3; ++k){
		
		temp_position = index-1+k;
		
		position[k] = temp_position;

		if(temp_position < 0){
			position[k] += numel-1;
		}else if(temp_position >= numel){
			position[k] -= numel-1;
		}
		
	}
	
	return 0;
}

int DensityProjector::get_coef(real_t *start, const int numel, const int shift){
	
	const double z1(sqrt(3)-2), tolerance(1e-15);
	const double lambda(6.);
	
	
	start[0] = get_first_causal(&start[0], numel, shift, tolerance);
	start[0] *= lambda;
	
	for(int k(1); k < numel; ++k){
		start[k*shift] *= lambda;
		start[k*shift] += z1 * start[ (k-1)*shift ];
	}
	
	// Last element of the line we are considering
	start[ (numel-1)*shift ] = get_last_anticausal(&start[0], numel, shift);
	
	for(int k(numel-2); k > -1; --k){
		start[k*shift] = z1*( start[ (k+1)*shift ] - start[ k*shift ] );
	}
	
	return 0;
}

real_t DensityProjector::get_first_causal(const real_t *start, const int numel, const int shift, const double tolerance){
	
	const double z1(sqrt(3)-2);
	const int k0( ceil( log(tolerance)/log(fabs(z1)) ) );
	double sum(start[0]);
	double zn(z1);
	
	const int num_sum = ( k0 < numel ? k0 : numel);
	
	for(int i(1); i < num_sum; ++i){
		sum += zn * start[i * shift];
		zn *= z1;
	}
		
	return sum;
}

real_t DensityProjector::get_last_anticausal(const real_t *start, const int numel, const int shift){
	
	const double z1(sqrt(3)-2);
	
	return - z1 / (1-z1*z1) * (start[ (numel-1) * shift ] + z1*start[ shift ]); 
	
}


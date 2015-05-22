#include "density_projector3d.hpp"

DensityProjector::DensityProjector(const int nx, const int ny, const int nz, double* data):
nx(nx), ny(ny), nz(nz), data(data){
	
	counter = new AtomicCounter();
	
	// std::string filename("data/cube128_1gaussian_1_1.bin");
	//
	// std::ifstream file;
	//
	// file.open(filename, std::ios::in | std::ios::binary);
	//
	//   	if (file.is_open()){
	//     	file.read(reinterpret_cast<char*>(&data[0]), nx * ny * nz * (sizeof data[0]) );
	//    	} else{
	//    		std::cout << "File " << filename <<" not opened" << std::endl;
	//    	}
	//

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

double DensityProjector::operator()(const madness::coord_3d& x) const{
	
	counter->increment();

	double 	returnvalue;
	double  interx, intery;			// Temporary variable to store the interpolated values along x and y axis
	int 	disp;					// Position

	std::vector<double> xweights(4);
	std::vector<double> yweights(4);
	std::vector<double> zweights(4);

	std::vector<int> xpositions(4);
	std::vector<int> ypositions(4);
	std::vector<int> zpositions(4);


	get_weights_and_position(Axis::X, x, xweights, xpositions);
	get_weights_and_position(Axis::Y, x, yweights, ypositions);
	get_weights_and_position(Axis::Z, x, zweights, zpositions);

	//printf("Weights %f\t %f\t %f\t %f\n", xweights[0], xweights[1], xweights[2], xweights[3]);

	returnvalue = 0.0;

	// Actually compute the value
	for(int k(0); k <= 3; ++k){

		intery = 0.0; // Interpolated value along x-axis

		for(int j(0); j <= 3; ++j){

			interx = 0.0; // Interpolated value along x-axis

			for(int i(0); i <= 3; ++i){
				disp = xpositions[i] + ypositions[j]*nx + zpositions[k]*nx*ny; // Position

				interx += xweights[i] * data[disp];
				// interx += xweights[i];
			}


			// if(std::abs(interx - 1.0) > 1e-14){
			// 	printf("\n interx is not 1, value is %f and error is %e\n", interx, std::abs(interx - 1.0));
			// 	// exit(0);
			// }
			
			
			
			intery += yweights[j] * interx;
			// intery should be one
		}
		returnvalue += zweights[k] * intery;
	}

	//printf("rho(%f, %f, %f) = %f\n",x[0], x[1], x[2], returnvalue);
	
	// if(std::abs(returnvalue - 1.0) > 1e-14){
	// 	printf("\nReturnvalue is not 1, value is %f and error is %e", returnvalue, std::abs(returnvalue - 1.0));
	// }
	//
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

int DensityProjector::get_weights_and_position(Axis axis, const madness::coord_3d& x, std::vector<double>& weights, std::vector<int>& position) const{
	
	int pos, index, numel;
	double alpha;
	
	numel = ( (axis == Axis::X) ? nx : ( (axis == Axis::Y) ? ny : nz) );
	pos = ( (axis == Axis::X) ? 0 : ( (axis == Axis::Y) ? 1 : 2) );
	
	index = (int) floor(x[pos] - 1);
	
	alpha = x[pos] - (index + 1); 	
	
	weights[0] = 1./6. * (1. - alpha) * (1. - alpha) * (1. - alpha);
	weights[1] = 2./3. - .5 * alpha * alpha * (2. - alpha);
	weights[2] = 2./3. - .5 * (1. - alpha) * (1. - alpha) * (1. + alpha);
	weights[3] = 1./6. * alpha * alpha * alpha;
	
	if (std::abs(weights[0] + weights[1] + weights[2] + weights[3] - 1.0) > 1e-13 ){
		printf("\nSum of weights is not equal to 1\n\n");
	}
	
	for(int k(0); k <= 3; ++k){
		
		position[k] = index - 1 + k;
		
		if(position[k] < 0){
			position[k] += numel-1; 
		}else if(position[k] >= numel){
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
	double sum(start[0]); // s(0)*z_1^0 = s(0)
	double zn(z1);
	
	const int num_sum = ( k0 < numel ? k0 : numel);
	
	for(int i(1); i < num_sum; ++i){
		sum += zn * start[i * shift];
		zn *= z1;
	}
	
	// if(k0 < numel){
	// 	for(int i(1); i < k0; ++i){
	// 		sum += zn * start[i * shift];
	// 		zn *= z1;
	// 	}
	// }
	// Add computation for "complete" loop
	
	
	return sum;
}

real_t DensityProjector::get_last_anticausal(const real_t *start, const int numel, const int shift){
	
	const double z1(sqrt(3)-2);
	
	return - z1 / (1-z1*z1) * (start[ (numel-1) * shift ] + z1*start[ shift ]); 
	
}


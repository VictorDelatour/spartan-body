#include <madness/mra/mra.h>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>


typedef double real_t;

struct AtomicCounter {
    std::atomic<int> value;

    void increment(){
        ++value;
    }

    void decrement(){
        --value;
    }

    int get(){
        return value.load();
    }
	
	void reset(){
		value = 0;
	}
};

class DensityProjector: public madness::FunctionFunctorInterface<double,3>{

public:
	
	DensityProjector(const int nx, const int ny, const int nz, double* data);
	~DensityProjector();
	
	double operator()(const madness::coord_3d& x) const;
	
	const int get_counter() const;
	void reset_counter();
	
	const double test_performance(const madness::coord_3d& x, const int& npoints) const;
	
	
	
protected:
	
private:
	
	enum class Axis { X, Y, Z };
	
	const int nx;
	const int ny;
	const int nz;
	
	real_t* data;
	
	// std::atomic<int> counter;
	AtomicCounter* counter;
	

	

	int inplace_filtering(Axis axis);
	
	int gaussian_filtering(Axis axis);
	
	int get_weights_and_position(Axis axis, const madness::coord_3d& x, std::vector<real_t>& weights, std::vector<int>& position) const;
	
	///
	/// @brief 	Transforms \a start into cubic B-spline interpolation coefficients
	///
	/// @param[in, out]	start	Pointer to the first element of the array to be processed
	/// @param[in]		numel	Number of elements of the array to be processed
	/// @param[in]		shift	Distance between two elements in the array (1 = x-axis, nx = y-axis, nx*ny = z-axis)
	///
	/// get_coef() can be applied on a row (shift = 1), on a column (shift = @a nx) or on a tower (shift = @a nx x @a ny)
	/// Independently of the shift chosen, the data is processed in place by applying a causal and an anti-causal filtering
	/// The causal filtering is initialised by a call to get_first_causal(), and the array is then processed in increasing 
	/// order, starting from start[0]. The anti-causal filtering is initialised by a call to get_last_anticausal(), 
	/// and the array is then processed in decreasing order, start from start[numel-1]
	///
	/// The algorithm for the causal and anti-causal filtering can be found in [Unser 99, p.26]:
	/// " Given the input signal values ${s(k)}_{k=0,\ldots, N-1}$ and defining $c^-(k)=\frac{c(k)}{6}$, 
	/// the right-hand-side factorization leads to the following recursive algorithm:
	/// $c^+(k) = s(k) + z_1 c^+(k-1)$			with $k = 1,\ldots,N-1$
	/// $c^-(k) = z_1 \left( c^-(k+1)-c^+(k) )$	with $k = N-2,\ldots,0$ "
	/// 
	/// Where $z_1 = \sqrt{3}-2$
	///
	int get_coef(real_t *start, const int numel, const int shift);
	
	///
	/// @brief Returns the initialisation of the causal filtering
	///
	/// @param[in]	start		Pointer to the first element of the array to be processed
	/// @param[in]	numel		Number of elements of the array to be processed
	/// @param[in]	shift 		Distance between two elements in the array (1 = x-axis, nx = y-axis, nx*ny = z-axis)
	/// @param[in]	tolerance	Tolerance used to compute the horizon of the causal initialization, default 1e-6
	///
	/// The initialisation for the causal filtering is given in [Unser 99, p.26], by
	/// $c^+(0) = \sum_{k=0}^{\infty}s(k)z_1^k$, where $z_1$ is the pole $\sqrt{3}-2$.
	/// In practice, the number of coefficients computer is $k_0$, with 
	/// $k_0 > \frac{ \log\varepsilon }{ \log \vert z_1\vert }$. 
	///
	/// Some authors [Ruijters and Thevenaz 2002] argue that $k_0 = 12$ is enough.
	///
	real_t get_first_causal(const real_t *start, const int numel, const int shift, const double tolerance);
	
	///
 	/// @brief Returns the initialisation of the anti-causal filtering
	/// 
	/// @param[in]	start		Pointer to the first element of the array to be processed
	/// @param[in]	numel		Number of elements of the array to be processed
	/// @param[in]	shift 		Distance between two elements in the array (1 = x-axis, nx = y-axis, nx*ny = z-axis)
	///
	/// The initialisation for the anti-causal filtering is given in [Unser 99, p.26], by
	/// c^-(N-1) = \frac{z_1}{1-z_1^2}\left( c^+(N-1) + z_1 c^+(N-2) )
	///
	real_t get_last_anticausal(const real_t *start, const int numel, const int shift);
	
};
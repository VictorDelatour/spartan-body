# spartan-body

Spartan-body uses the [MADNESS](https://github.com/m-a-d-n-e-s-s/madness/) numerical environment to solve the Poisson equation for the gravitational potential and was implemented as part of an internship. The goal of said internship was to test the use of MADNESS as a viable alternative Poisson solver and compare its accuracy and performance with a particle-mesh method based Poisson solver using the Fastest Fourier Transform in the West (FFTW) library.

## Usage

### Dependencies

To use the the code, you'll need the [MADNESS](https://github.com/m-a-d-n-e-s-s/madness/) library on your computer. Download it and build it, everything is explained in details [here](https://github.com/m-a-d-n-e-s-s/madness/wiki#Building_MADNESS "Link to madness build directives"). When `make` is finished, you can find an example makefile for applications using madness.

### Building and running the executable

Using the provided makefile, simply `make` the spartan-body executable using GNU compilers for Fortran and C++.

To run the executable for a gridsize *nx*
```
arpun -n #NUM_OF_PROCS -cc none ./spartan nx
```
failing to remove thread pinning ('-cc none') reduces drastically the performance. 

#### Data

The particles used originate from a RAMSES simulation for astrophysics, in order to provide a realistic initialization, and are expected to be stored in the *output_00003* folder. At the moment, the only allowed size *nx* allowed are 128, 256 and 512. Based on the variable *nx*, the code will read an `info_0000*.txt` file containing the number of CPUs used to produce the particles, which corresponds to the number of files used to store the particles. Each `part_0000*.out0000#CPU` is expected to contain the number of particles contained in the file, as well as their initial position, velocity and mass.


## References

Here's a list of all the references and documentation that were used during this project.  

### madness and multiwavelets

Mathematical background for the madness environment

**Alpert, Beylkin, et al.** ["Adaptive solution of partial differential equations in multiwavelet bases."](http://math.nist.gov/~BAlpert/mwpde.pdf "Adaptive solution of partial differential equations in multiwavelet bases") Journal of Computational Physics 182.1 (2002): 149-190.

**Harrison, Robert J., et al.** ["Multiresolution quantum chemistry: Basic theory and initial applications."](http://amath.colorado.edu/faculty/beylkin/papers/H-F-Y-G-B-2004.pdf "Multiresolution quantum chemistry: Basic theory and initial applications") The Journal of chemical physics 121.23 (2004): 11587-11598.

**Thornton, W. Scott, Nicholas Vence, and Robert Harrison.** ["Introducing the MADNESS numerical framework for petascale computing."](https://scholar.google.ch/scholar?q=Introducing+the+MADNESS+numerical+framework+for+petascale+computing&btnG=&hl=en&as_sdt=0%2C5 "Introducing the MADNESS numerical framework for petascale computing.") Proceedings of the Cray Users Group (2009).

### Cubic B-spline and interpolation

Mathematical theory behind cubic B-spline and implementation details of interpolation

**Ruijters, Daniel, and Philippe Thévenaz.** ["GPU prefilter for accurate cubic B-spline interpolation."](http://bigwww.epfl.ch/publications/ruijters1201.pdf "GPU prefilter for accurate cubic B-spline interpolation") The Computer Journal (2010): bxq086.

**Thévenaz, Philippe, Thierry Blu, and Michael Unser.** ["Interpolation revisited."](http://infoscience.epfl.ch/record/63069/files/thevenaz0002.pdf "Interpolation revisited") Medical Imaging, IEEE Transactions on 19.7 (2000): 739-758.

**Ruijters, Daniel, Bart M. ter Haar Romeny, and Paul Suetens.** ["Efficient GPU-based texture interpolation using uniform B-splines."](http://www.mate.tue.nl/mate/pdfs/10318.pdf "Efficient GPU-based texture interpolation using uniform B-splines") Journal of Graphics, GPU, and Game Tools 13.4 (2008): 61-69.

**Unser, Michael.** ["Splines: A perfect fit for signal and image processing."](http://infoscience.epfl.ch/record/63064/files/unser9902.pdf "Splines: A perfect fit for signal and image processing.") Signal Processing Magazine, IEEE 16.6 (1999): 22-38.

**Unser, Michael, Akram Aldroubi, and Murray Eden.** ["B-spline signal processing. I. Theory."](http://bigwww.epfl.ch/publications/unser9301.pdf "B-spline signal processing. I. Theory") Signal Processing, IEEE Transactions on 41.2 (1993): 821-833.

**Unser, Michael, Akram Aldroubi, and Murray Eden.** ["B-spline signal processing. II. Efficiency design and applications."](http://users.fmrib.ox.ac.uk/~jesper/papers/future_readgroups/unser9302.pdf "B-spline signal processing. II. Efficiency design and applications") Signal Processing, IEEE Transactions on 41.2 (1993): 834-848.

**Unser, Michael, Akram Aldroubi, and Murray Eden.** ["Fast B-spline transforms for continuous image representation and interpolation."](http://bigwww.epfl.ch/publications/unser9102.pdf "Fast B-spline transforms for continuous image representation and interpolation") IEEE Transactions on Pattern Analysis & Machine Intelligence 3 (1991): 277-285.

### Scatered data interpolation

Algorithms and implementation of scatered data interpolation. This hasn't been implemented in the spartan-body code, but could be a solution to cirumvent the projection of the particles on the uniform cubic mesh and thus improve the accuracy of the method. A possible solution was the use of RBF or Shepard's algorithm

**Anjyo, Ken, J. P. Lewis, and Frédéric Pighin**. ["Scattered data interpolation for computer graphics."](http://scribblethink.org/Courses/ScatteredInterpolation/scatteredinterpcoursenotes.pdf) ACM SIGGRAPH 2014 Courses. ACM, 2014.

**Amidror, Isaac.** ["Scattered data interpolation methods for electronic imaging systems: a survey."](http://infoscience.epfl.ch/record/99883/files/sdimfeisas.pdf) Journal of electronic imaging 11.2 (2002): 157-176.

**Awanou, Gerard, Ming-Jun Lai, and Paul Wenston.** ["The multivariate spline method for scattered data fitting and numerical solutions of partial differential equations."](http://homepages.math.uic.edu/~awanou/multi.pdf) Wavelets and splines: Athens (2005): 24-74.

**Bertram, Martin, Xavier Tricoche, and Hans Hagen.** ["Adaptive smooth scattered-data approximation for large-scale terrain visualization."](https://www.cs.purdue.edu/homes/xmt/papers/terrain.pdf) VisSym. 2003.

**Feng, Renzhong, and Yanan Zhang.** ["Piecewise Bivariate Hermite Interpolations for Large Sets of Scattered Data."](http://www.emis.de/journals/HOA/JAM/Volume2013/239703.pdf) Journal of Applied Mathematics 2013 (2013).

**Lazzaro, Damiana.** "A parallel multivariate interpolation algorithm with radial basis functions." International journal of computer mathematics 80.7 (2003): 907-919.

**Lazzaro, Damiana, and Laura B. Montefusco.** ["Radial basis functions for the multivariate interpolation of large scattered data sets."](http://www.sciencedirect.com/science/article/pii/S037704270100485X) Journal of Computational and Applied Mathematics 140.1 (2002): 521-536.

**Lee, Seungyong, George Wolberg, and Sung Yong Shin.** ["Scattered data interpolation with multilevel B-splines."](http://www.researchgate.net/profile/George_Wolberg/publication/3410822_Scattered_data_interpolation_with_multilevel_B-splines/links/00b49518719ac9f08a000000.pdf) Visualization and Computer Graphics, IEEE Transactions on 3.3 (1997): 228-244.

**Renka, Robert J.** "Multivariate interpolation of large sets of scattered data." ACM Transactions on Mathematical Software (TOMS) 14.2 (1988): 139-148.

**Renka, Robert J.** "Algorithm 660: QSHEP2D: Quadratic Shepard method for bivariate interpolation of scattered data." ACM Transactions on Mathematical Software (TOMS) 14.2 (1988): 149-150.

**Thacker, William I., et al.** ["Algorithm 905: SHEPPACK: Modified Shepard algorithm for interpolation of scattered multivariate data."](http://s3.amazonaws.com/researchcompendia_prod/articles/5aba70101f0bbfbfc8733cbceee19109-2013-12-23-01-48-28/a34-thacker.pdf) ACM Transactions on Mathematical Software (TOMS) 37.3 (2010): 34.


















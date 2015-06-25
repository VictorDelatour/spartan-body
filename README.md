# spartan-body

Spartan-body is an n-body solver that uses the [MADNESS](https://github.com/m-a-d-n-e-s-s/madness/) numerical environment and was written as part of an internship. The goal of said internship was to test the feasibility of an n-body solver based on MADNESS and compare its accuracy and performance with a particle-mesh method based n-body solver using the Fastest Fourier Transform in the West library.

## Usage

## References

Here's a list of all the references that were used during the internship

### madness and multiwavelets

**Alpert, Beylkin, et al.** ["Adaptive solution of partial differential equations in multiwavelet bases."](http://math.nist.gov/~BAlpert/mwpde.pdf "Adaptive solution of partial differential equations in multiwavelet bases") Journal of Computational Physics 182.1 (2002): 149-190.

**Harrison, Robert J., et al.** ["Multiresolution quantum chemistry: Basic theory and initial applications."](http://amath.colorado.edu/faculty/beylkin/papers/H-F-Y-G-B-2004.pdf "Multiresolution quantum chemistry: Basic theory and initial applications") The Journal of chemical physics 121.23 (2004): 11587-11598.

**Thornton, W. Scott, Nicholas Vence, and Robert Harrison.** ["Introducing the MADNESS numerical framework for petascale computing."](https://scholar.google.ch/scholar?q=Introducing+the+MADNESS+numerical+framework+for+petascale+computing&btnG=&hl=en&as_sdt=0%2C5 "Introducing the MADNESS numerical framework for petascale computing.") Proceedings of the Cray Users Group (2009).

### Cubic B-spline and interpolation


**Ruijters, Daniel, and Philippe Thévenaz.** ["GPU prefilter for accurate cubic B-spline interpolation."](http://bigwww.epfl.ch/publications/ruijters1201.pdf "GPU prefilter for accurate cubic B-spline interpolation") The Computer Journal (2010): bxq086.

**Thévenaz, Philippe, Thierry Blu, and Michael Unser.** ["Interpolation revisited."](http://infoscience.epfl.ch/record/63069/files/thevenaz0002.pdf "Interpolation revisited") Medical Imaging, IEEE Transactions on 19.7 (2000): 739-758.

**Ruijters, Daniel, Bart M. ter Haar Romeny, and Paul Suetens.** ["Efficient GPU-based texture interpolation using uniform B-splines."](http://www.mate.tue.nl/mate/pdfs/10318.pdf "Efficient GPU-based texture interpolation using uniform B-splines") Journal of Graphics, GPU, and Game Tools 13.4 (2008): 61-69.

**Unser, Michael.** ["Splines: A perfect fit for signal and image processing."](http://infoscience.epfl.ch/record/63064/files/unser9902.pdf "Splines: A perfect fit for signal and image processing.") Signal Processing Magazine, IEEE 16.6 (1999): 22-38.

**Unser, Michael, Akram Aldroubi, and Murray Eden.** ["B-spline signal processing. I. Theory."](http://bigwww.epfl.ch/publications/unser9301.pdf "B-spline signal processing. I. Theory") Signal Processing, IEEE Transactions on 41.2 (1993): 821-833.

**Unser, Michael, Akram Aldroubi, and Murray Eden.** ["B-spline signal processing. II. Efficiency design and applications."](http://users.fmrib.ox.ac.uk/~jesper/papers/future_readgroups/unser9302.pdf "B-spline signal processing. II. Efficiency design and applications") Signal Processing, IEEE Transactions on 41.2 (1993): 834-848.

**Unser, Michael, Akram Aldroubi, and Murray Eden.** ["Fast B-spline transforms for continuous image representation and interpolation."](http://bigwww.epfl.ch/publications/unser9102.pdf "Fast B-spline transforms for continuous image representation and interpolation") IEEE Transactions on Pattern Analysis & Machine Intelligence 3 (1991): 277-285.

### Scatered data interpolation

**Amidror, Isaac.** ["Scattered data interpolation methods for electronic imaging systems: a survey."](http://infoscience.epfl.ch/record/99883/files/sdimfeisas.pdf) Journal of electronic imaging 11.2 (2002): 157-176.

**Awanou, Gerard, Ming-Jun Lai, and Paul Wenston.** ["The multivariate spline method for scattered data fitting and numerical solutions of partial differential equations."](http://homepages.math.uic.edu/~awanou/multi.pdf) Wavelets and splines: Athens (2005): 24-74.

**Lee, Seungyong, George Wolberg, and Sung Yong Shin.** ["Scattered data interpolation with multilevel B-splines."](http://www.researchgate.net/profile/George_Wolberg/publication/3410822_Scattered_data_interpolation_with_multilevel_B-splines/links/00b49518719ac9f08a000000.pdf) Visualization and Computer Graphics, IEEE Transactions on 3.3 (1997): 228-244.

**Renka, Robert J.** "Multivariate interpolation of large sets of scattered data." ACM Transactions on Mathematical Software (TOMS) 14.2 (1988): 139-148.

**Thacker, William I., et al.** ["Algorithm 905: SHEPPACK: Modified Shepard algorithm for interpolation of scattered multivariate data."](http://s3.amazonaws.com/researchcompendia_prod/articles/5aba70101f0bbfbfc8733cbceee19109-2013-12-23-01-48-28/a34-thacker.pdf) ACM Transactions on Mathematical Software (TOMS) 37.3 (2010): 34.

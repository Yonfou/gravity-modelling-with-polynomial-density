# gravity-modelling-with-polynomial-density

## 3D Large-scale Gravitational Modelling Based on Tesseroids with a Polynomial Density of Arbitrary Degree in depth

by Fang Ouyang and Longwei Chen

The algorithm is aimed at achieving the 3D large-scale gravitational forward modelling in an efficient and accurate way. It can deal with a polynomial density up to an arbitrary order in depth. In the algorithm, the source region is divided into a number of tesseroids, and the density in each tesseroid is assumed to be a polynomial function of arbitrary degree. To guarantee the computational accuracy and accelerate the calculation, two key points are involved: (1) the volume Newton's integral is decomposed into a one-dimensional integral with a polynomial density in the radial direction, which is solved analytically in a recursive way, and a surface integral over the horizontal directions evaluated by the Gaussian Legendre quadrature (GLQ) combined with a 2D adaptive discretization; (2) a fast and flexible discrete convolution algorithm based on 1D fast Fourier Transform (FFT) and a general Toepritz form of weight coefficient matrices is adopted in the longitudinal dimension to speed up the computation of the cumulative contributions from all tesseroids.

## How to run the codes?

We provide a shell model for reference.

Step 1:  Establish a density model and generalize the model files required in 'Para.txt'. 
            For example, the Matlab file 'shell\cook_shell.m'  helps one to established a shell model.
Step 2: Modify the file path of 'Para.txt' in the source code 'ReadIn.f90' according to your own needs.
Step 3: Run the codes.


## License

All source code is made available under a BSD 3-clause license. You can freely use and modify the code, without warranty, so long as you provide attribution to the authors. See LICENSE.md for the full license text.

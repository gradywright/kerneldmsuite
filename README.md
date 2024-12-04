# Kernel DM Suite
A MATLAB suite of methods for constructing differentiation matrices (DMs) on point clouds using kernel methods.  

Techniques supported so far:

* Sphere (see [[1]](#HRW24)) for details)
    * Radial Basis Function Finite Difference (RBF-FD)
    * Local Lagrange 
    * Global or Pseudospectral

All of the techniques support various kernels and the inclusion of polynomials in the approximation space.

# Installation and requirements

The Kernel DM Suite is compatible with MATLAB R2018a and later.  It requires that the MATLAB statistics toolbox be installed because it makes use of that toolboxes KD-tree.

To install, clone the directory with Git:
```
git clone https://github.com/gradywright/kerneldmsuite.git
```
You will then need to add the `kerneldmsuite` directory to the MATLAB path:
```
addpath(kdmroot), savepath
```
where `kdmroot` is the path to the kerneldmsuite directory

# Getting started

The following example shows how to use the code to produce different DMs for the Laplace-Beltrami operator on the sphere.  This example uses the Hammersley point set on the sphere and the restricted surface spline kernel of order 2 ($`\phi(r) = r^4 \log(r)`$) with degree 2 spherical harmonics precision.
```matlab
% Point cloud - Hammersley points
N = 4096;
x = kdmutil.sphHammersleyNodes(N);
% Kernel: r^4*log(r)
rbf = kernel.phseven(2);
% Stencil size
n = 41;
% Spherical harmonic degree 
deg = 2;
% DM using RBF-FD method
L = kdmsphere.rbffdlap(x,rbf,n,deg);
% DM using Local Lagrange method
L = kdmsphere.loclaglap(x,rbf,n,deg);
% DM using Global method
L = kdmsphere.globrbflap(x,rbf,deg);
```
More examples can be found in [examples](examples/) folder.

To create DMs based on different point clouds on the sphere, use the [spherepts](https://github.com/gradywright/spherepts) package.

# References:

<a name="HRW24">[1]</a> T. Hangelbroek, C. Rieger, and G. B. Wright. Spectral stability and perturbation results for kernel differentiation matrices on the sphere. [arXiv:2311.06982](https://arxiv.org/abs/2311.06982)


# Acknowledgements 
This software development was partially supported by National Science Foundation grant 2309712.


https://github.com/gradywright/kerneldmsuite#references





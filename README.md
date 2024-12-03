# Kernel DM Suite
A MATLAB suite of methods for constructing differentiation matrices on point clouds computed using kernel methods.  

Techniques supported so far:

* Sphere (see [1] for details)
+ Radial Basis Function Finite Difference (RBF-FD)
+ Local Lagrange 
+ Global or Pseudospectral

All of the techniques support various kernels and the inclusion of polynomials in the approximation space.

The following example shows how to use the code to produce approximations to the Laplace-Beltrami operator on the sphere using the restricted surface spline kernel of order 2 with degree 2 spherical harmonics precision.
```
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
L = kdmsphere.globlap(x,rbf,deg);
```
More examples can be found in [examples.m](https://raw.github.com/gradywright/kerneldmsuite/master/examples.m).

For differentiation matrices based on different point clouds on the sphere, use the [spherepts](https://github.com/gradywright/spherepts) package.

## References:

[1] T. Hangelbroek, C. Rieger, and G. B. Wright. Spectral stability and perturbation results for kernel differentiation matrices on the sphere. [arXiv:2311.06982](https://arxiv.org/abs/2311.06982)


## Acknowledgements 
This software development was partially supported by National Science Foundation grant 2309712.








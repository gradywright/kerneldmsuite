%% Differentiation matrices for the Laplace-Beltrami operator on the sphere
% The code below gives several examples for approximating the spherical
% Laplacian of the a function sampled a a point cloud using kernel
% differentiation matrices.

%
% Set-up parameters for constructing the differentiation matrices.
%

% Use Hammersley points.  More options are available from the spherepts package
N = 4096;
x = kdmutil.sphHammersleyNodes(N);

% Use the restricted polyharmonic spline (also called the surface spline) of
% order 2 (i.e., phi(r) = r^4*log(r)
order = 2;
rbf = kernel.phseven(order);

% Append spherical harmonics of degree 2
deg = 2;

% Formula for determining the stencil size for RBF-FD and Local Lagrange methods
stencilSize = @(K) ceil(K.^2*log(N)^2/7);

% Use a Gaussain to test the resulting approximations of the Laplace-Beltrami
xc = [1 0 0];  % Center of the Gaussian
r2 = (xc(1,1)-x(:,1)).^2 + (xc(1,2)-x(:,2)).^2 + (xc(1,3)-x(:,3)).^2;
sig2 = 1;
f = exp(-sig2 * r2);
exactlap = sig2*exp(-sig2*r2).*(-4 + r2.*(2 - sig2*(-4 + r2)));


%% Differentiation matrix using the RBF-FD method
% Stencil size
n = stencilSize(2);
fprintf('Construcing RBF-FD DM with stencil size = %d...\n',n)

tic
L = kdmsphere.rbffdlap(x,rbf,n,deg);
etime = toc;

fprintf('Finished in %1.3e s\n',etime);

% Apply to f and compute the error
lapf = L*f;
err = norm(lapf - exactlap,inf);

fprintf('Max-norm error = %.4e\n\n',err)

%% Differentiation matrix using the Local-Lagrange method
% Stencil size
n = stencilSize(7);
fprintf('Construcing Local Lagrange DM with stencil size = %d...\n',n)

tic
L = kdmsphere.loclaglap(x,rbf,n,deg);
etime = toc;

fprintf('Finished in %1.3e s\n',etime);

% Apply to f and compute the error
lapf = L*f;
err = norm(lapf - exactlap,inf);

fprintf('Max-norm error = %.4e\n\n',err)

%% Differentiation matrix using the Global (or pseudospectral) RBF method
fprintf('Construcing Global RBF DM...\n')

tic
L = kdmsphere.globrbflap(x,rbf,deg);
etime = toc;

fprintf('Finished in %1.3e s\n',etime);

% Apply to f and compute the error
lapf = L*f;
err = norm(lapf - exactlap,inf);

fprintf('Max-norm error = %.4e\n\n',err)


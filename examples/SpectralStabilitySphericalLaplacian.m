%% Spectral stability of kernel differentiation matrices (DMs) on the sphere
% The code below gives several examples from the paper:
%    T. Hangelbroek, C. Rieger, and G. B. Wright. Spectral stability and
%    perturbation results for kernel differentiation matrices on the sphere.
%    arXiv:2311.06982 
% This example file requires the spherepts package
% (https://github.com/gradywright/spherepts)

%
% Set-up parameters for constructing the differentiation matrices (DMs) for the
% Laplace-Beltrami operator on the sphere
%

% Use minimum energy points from the spherepts package
N = 2048;
x = getMinEnergyNodes(N);

% Use the restricted polyharmonic spline (also called the surface spline) of
% order 2 (i.e., phi(r) = r^4*log(r)
order = 2;
rbf = kernel.phseven(order);

% Append spherical harmonics of degree 2
deg = 2;

% Formula for determining the stencil size for RBF-FD and Local Lagrange methods
stencilSize = @(K) ceil(K.^2*log(N)^2/7);

% Plotting variables
FS = 'FontSize';
fs = 14;
INTERP = 'Interpreter';
interp = 'latex';
MS = 'MarkerSize';
ms = 6;

%% Spectra of the DM using the Global (or pseudospectral) RBF method
fprintf('Construcing Global RBF DM...\n')

tic
L = kdmsphere.globrbflap(x,rbf,deg);
etime = toc;

fprintf('Finished in %1.3e s\n',etime);

% Take the negative of L to correspond to the definition used in the paper.
L = -L;

% Compute the full spectrum
evglob = eig(L);

% Plot the spectrum
plot(real(evglob),imag(evglob),'rx',MS,ms)
xlabel('Re$(\mu)$',INTERP,interp,FS,fs)
ylabel('Im$(\mu)$',INTERP,interp,FS,fs)
title('Spectrum of the global RBF DM');

%% Spectra of the DM using the RBF-FD method
% Stencil size
n = stencilSize(7);
fprintf('Construcing RBF-FD DM with stencil size = %d...\n',n)

tic
L = kdmsphere.rbffdlap(x,rbf,n,deg);
etime = toc;

fprintf('Finished in %1.3e s\n',etime);

% Take the negative of L to correspond to the definition used in the paper.
L = -L;

% Compute the full spectrum
evrbffd = eig(full(L));

% Plot the spectrum and compare to global method
plot(real(evglob),imag(evglob),'rx',MS,ms), hold on
plot(real(evrbffd),imag(evrbffd),'bs',MS,ms)
xlabel('Re$(\mu)$',INTERP,interp,FS,fs)
ylabel('Im$(\mu)$',INTERP,interp,FS,fs)
title('Spectrum of the RBF-FD DM');
legend('Global RBF','RBF-FD',FS,fs)
hold off

%% Differentiation matrix using the Local-Lagrange method
% Stencil size
n = stencilSize(7);
fprintf('Construcing Local Lagrange DM with stencil size = %d...\n',n)

tic
L = kdmsphere.loclaglap(x,rbf,n,deg);
etime = toc;

fprintf('Finished in %1.3e s\n',etime);

% Take the negative of L to correspond to the definition used in the paper.
L = -L;

% Compute the full spectrum
evloclag = eig(full(L));

% Plot the spectrum and compare to global method
plot(real(evglob),imag(evglob),'rx',MS,ms), hold on
plot(real(evloclag),imag(evloclag),'bs',MS,ms)
xlabel('Re$(\mu)$',INTERP,interp,FS,fs)
ylabel('Im$(\mu)$',INTERP,interp,FS,fs)
title('Spectrum of the RBF-FD DM');
legend('Global RBF','Local Lagrange',FS,fs)
hold off
function Lap = globrbflap(nodes,rbf,deg)
% GLOBLAP Differentiation matrix for the Laplacian on the unit sphere using the
% global RBF method, or pseudospectral method.
%
% LapDM = globlap(nodes,rbf,deg) returns the differentiation
% matrix (DM) for the Laplace-Beltrami operator on the sphere using the global
% method, or pseudospectral method, which uses all the points on the sphere
% in forming the approximation. Parameters are as follows:
% nodes:   Locations on the unit sphere where the Laplace-Beltrami operator will
%          be approximated
% rbf:     RBF kernel object to use for the approximations
% deg:     Spherical harmonic degree of precision of the formulas, where -1
%          means no precision
%
% see also RBFFDLAP and LOCLAGLAP

% Copyright 2024 by Grady B. Wright

% TODO: Add better error checking

if ~isa(rbf,'kernel.rbf')
    error('KDMSUITE:DMSPHERE:globrbflap:rbf','The second input must be an RBF object.  For example, p = kernel.phsodd(2) or p = kernel.phseven(2).')
end

deg = round(deg);
if deg < -1
    error('KDMSUITE:DMSPHERE:globrbflap:degree','The spherical harmonic degree of precsion must be an integer >= -1.')
end

% Dimension of the space of spherical harmonics that will be used.
L = (deg+1)^2;

[N,d] = size(nodes);

if d ~= 3
    error('KDMSUITE:DMSPHERE:globrbflap:dimension','The nodes array should be of size N-by-3, where each of the N rows corresponds to a point in 3-dimensional space on the sphere.')
end

% Compute the eigenvalues for spherical harmonics
evlap = zeros(1,L);
degSph = zeros(L,1);
cnt = 1;
ZM = zeros(L);
for l=0:deg
    for k=-l:l
        evlap(cnt) = -l*(l+1);
        degSph(cnt) = l;
        cnt = cnt + 1;
    end
end

% Turn off the annoying warnings
warnstate = warning('query','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:nearlySingularMatrix');

[xk,xj] = meshgrid(nodes(:,1));
[yk,yj] = meshgrid(nodes(:,2));
[zk,zj] = meshgrid(nodes(:,3));
rd2 = ((xj-xk).^2 + (yj-yk).^2 + (zj-zk).^2);
rd = sqrt(rd2);

% Local kernel matrix
Phi = rbf.phi(rd);
% Laplacian of the kernel shifted to the stencil center and evaluated at
% all the stencil points.
Blap = rd2.*(4-rd2).*rbf.zeta(rd)/4 + (2-rd2).*rbf.eta(rd);

if deg > -1
    P = kdmutil.sphHarmBasis(nodes,deg);
    A = [[Phi P];[P.' ZM]];
    Plap = P.*evlap;
    Lap = [Blap Plap]/A;
    Lap = Lap(1:N,1:N);
else
    Lap = Blap/Phi;
end

% Return the warning state
warning(warnstate);

end

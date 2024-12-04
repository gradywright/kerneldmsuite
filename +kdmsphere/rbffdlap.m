function Lap = rbffdlap(nodes,rbf,nd,deg)
% RBFFDLAP Differentiation matrix for the Laplacian on the unit sphere using the
% RBF-FD method
%
% LapDM = rbffddm(nodes,rbf,nd,deg) returns the differentiation
% matrix (DM) for the Laplace-Beltrami operator on the sphere using the RBF-FD
% method. Parameters are as follows:
% nodes:   Locations on the unit sphere where the Laplace-Beltrami operator will
%          be approximated
% rbf:     RBF kernel object to use for the approximations
% nd:      Size of the stencil to use
% deg:     Spherical harmonic degree of precision of the formulas, where -1
%          means no precision
%
% see also LOCLAGLAP and GLOBRBFLAP

% Copyright 2024 by Grady B. Wright

% TODO: Add better error checking

if ~isa(rbf,'kernel.rbf')
    error('KDMSUITE:DMSPHERE:rbffdlap:rbf','The second input must be an RBF object.  For example, p = kernel.phsodd(2) or p = kernel.phseven(2).')
end

deg = round(deg);
if deg < -1
    error('KDMSUITE:DMSPHERE:rbffdlap:degree','The spherical harmonic degree of precsion must be an integer >= -1.')
end

% Dimension of the space of spherical harmonics that will be used.
L = (deg+1)^2;

nd = round(nd);
if nd <= L
    error('KDMSUITE:DMSPHERE:rbffdlap:stencilSize','The stencil size needs to be larger than the dimension of the space of spherical harmonics of the given polyDeg, which is L=%d.  A common choice would be nd=%d.',L,ceil(1.5*L));
end

[N,d] = size(nodes);

if d ~= 3
    error('KDMSUITE:DMSPHERE:rbffdlap:dimension','The nodes array should be of size N-by-3, where each of the N rows corresponds to a point in 3-dimensional space on the sphere.')
end

% Set up the arrays to contain the weights for each stencil
row_index = repmat(1:N,[nd 1]);
col_index = row_index;
wghts_lap = row_index;

if ~license('test','statistics_toolbox')
    error('KDMSUITE:DMSPHERE:rbffdlap:kdtree','This code requires the statistics toolbox to build and search a KD-tree of the points.')
end
treeroot = createns(nodes);
[idx,~] = knnsearch(treeroot,nodes,'k',nd);

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

% Initialize code to run in parallel
if isempty(gcp('nocreate'))
    parpool('Threads');
end

% Turn off the annoying warnings
warnstate = warning('query','MATLAB:nearlySingularMatrix');
parfevalOnAll(@warning,0,'off','MATLAB:nearlySingularMatrix');

parfor i=1:N
    %
    % Extract out the points for the ith node.
    %
    j = idx(i,:);
    x = nodes(j,:);

    % This rotation stuff is probably not needed, but it does keep the
    % evaluation of the spherical harmonics consistent.

    % % Convert the current point to lattitude longitude
    % [lam,th] = cart2sph(x(1,1),x(1,2),x(1,3));
    % 
    % % Rotate the point to the north pole.
    % plam = lam;
    % pth = pi-th;
    % D = [[cos(plam) sin(plam) 0];[-sin(plam) cos(plam) 0];[0 0 1]];
    % C = [[sin(pth) 0 cos(pth)];[0 1 0];[-cos(pth) 0 sin(pth)]];    
    % % Rotate the center and all the other points.
    % x = (C*(D*x.')).';
    
    nd = length(x);
    [xk,xj] = meshgrid(x(:,1));
    [yk,yj] = meshgrid(x(:,2));
    [zk,zj] = meshgrid(x(:,3));
    rd2 = ((xj-xk).^2 + (yj-yk).^2 + (zj-zk).^2);
    rd = sqrt(rd2);
    
    % Local kernel matrix
    Phi = rbf.phi(rd);
    % Laplacian of the kernel shifted to the stencil center and evaluated at
    % all the stencil points.
    Blap = rd2(:,1).*(4-rd2(:,1)).*rbf.zeta(rd(:,1))/4 + (2-rd2(:,1)).*rbf.eta(rd(:,1));

    if deg > -1
        P = kdmutil.sphHarmBasis(x,deg);
        A = [[Phi P];[P.' ZM]];
        Plap = evlap.*P(1,:);
        Lap = A\[Blap;Plap.'];
    else
        Lap = Phi\Blap;
    end
    
    % Store the weights and the rows and columns that they correspond to.
    wghts_lap(:,i) = Lap(1:nd,1);
    row_index(:,i) = i*ones(1,nd);
    col_index(:,i) = j;
end
% Build the sparse differentiation matrix from the weights
Lap = sparse(row_index(:),col_index(:),wghts_lap(:),N,N);

% Return the warning state
warning(warnstate);

end

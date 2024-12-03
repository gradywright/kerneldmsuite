function pass = test_loclaglapdm( )

% Base tolerance for the tests
tol = 5e-4;

% Node set to test
N = 400;
x = kdmutil.sphHammersleyNodes(N);
N = size(x,1);

% Test that the method gives a DM of the correct size
stencilSize = 61;
deg = 2;
order = 2;
p = kernel.phseven(order);
L = kdmsphere.loclaglap(x,p,stencilSize,deg);

[NN,MM] = size(L);

j = 1;
pass(j) = (NN == N) && (MM == N);
j = j+1;

% Check that the matrix has the correct number of non-zeros
sparsity = nnz(L);
pass(j) = sparsity == N*stencilSize;
j = j+1;

stencilSize = ceil(0.95*N);
deg = 2;
order = 2;
p = kernel.phseven(order);
L = kdmsphere.loclaglap(x,p,stencilSize,deg);
sparsity = nnz(L);
pass(j) = sparsity == N*stencilSize;
j = j+1;

% Test that the method approximates spherical harmonics of degree 2 for large
% stencil sizes
egvls = [0 -2*ones(1,3) -6*ones(1,5)];
f = kdmutil.sphHarmBasis(x,2);
err = L*f - egvls.*f;
pass(j:j+8) = max(abs(err)) < tol;
j = j+9;

% Check that the method is exact for different RBFs.
rid = ceil(N/3);
xc = x(rid,:);
rd2 = (x(:,1)-xc(1)).^2 + (x(:,2)-xc(2)).^2 + (x(:,3)-xc(3)).^2;
rd = sqrt(rd2);
f = p.phi(rd(:,1));
lapf = rd2(:,1).*(4-rd2(:,1)).*p.zeta(rd(:,1))/4 + (2-rd2(:,1)).*p.eta(rd(:,1));
err = L(rid,:)*f - lapf(rid);

pass(j) = max(abs(err)) < 10*tol;
j = j+1;

% Check that the correct errors are thrown
try
    L = kdmsphere.loclaglap(x,p,stencilSize,-2);
    pass(j) = false;
catch ME
    pass(j) = strcmp(ME.identifier, 'KDMSUITE:DMSPHERE:loclaglap:degree');
end
j = j+1;

try
    L = kdmsphere.loclaglap(x,p,1,2);
    pass(j) = false;
catch ME
    pass(j) = strcmp(ME.identifier, 'KDMSUITE:DMSPHERE:loclaglap:stencilSize');
end
j = j+1;

try
    L = kdmsphere.loclaglap(x(:,1:2),p,stencilSize,2);
    pass(j) = false;
catch ME
    pass(j) = strcmp(ME.identifier, 'KDMSUITE:DMSPHERE:loclaglap:dimension');
end
j = j+1;

try
    L = kdmsphere.loclaglap(x,1,stencilSize,-1);
    pass(j) = false;
catch ME
    pass(j) = strcmp(ME.identifier, 'KDMSUITE:DMSPHERE:loclaglap:rbf');
end
j = j+1;

end


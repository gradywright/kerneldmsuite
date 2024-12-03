function P = sphHarmBasis(X,l)
% SPHHARMBASIS Evaluates each term in a standard basis of spherical harmonics 
% of on the unit sphere at the points in X.
%
% P = sphHarmBasis(X,l) evaluates L = (l+1)^2 spherical harmonic basis functions
% at the n values in X on the unit sphere.  P is an array of dimension n-by-L.
%
% Note that spherical harmonics up to degree l = 6 are supported.

% Copyright 2024 by Grady B. Wright

L = (l+1)^2;
P = zeros(size(X,1),L);
x = X(:,1);
y = X(:,2);
z = X(:,3);

% Assume a unit sphere
r2 = 1;

j = 1;
for k=0:l
    switch k
        case 0
            P(:,j) = 1+0*x;
            j = j + 1;
        case 1
            P(:,j) = y;
            j = j+1;
            P(:,j) = z;
            j = j+1;
            P(:,j) = x;
            j = j+1;
        case 2
            P(:,j) = x.*y;
            j = j + 1;
            P(:,j) = y.*z;
            j = j + 1;
            P(:,j) = 3*z.^2 - r2;
            j = j + 1;
            P(:,j) = x.*z;
            j = j + 1;
            P(:,j) = x.^2 - y.^2;
            j = j + 1;
        case 3
            P(:,j) = y.*(-3*x.^2 + y.^2);
            j = j + 1;
            P(:,j) = x.*y.*z;
            j = j + 1;
            P(:,j) = y.*(r2 - 5*z.^2);
            j = j + 1;
            P(:,j) = z.*(-3*r2 + 5*z.^2);
            j = j + 1;
            P(:,j) = x.*(r2 - 5*z.^2);
            j = j + 1;
            P(:,j) = z.*(x.^2 - y.^2);
            j = j + 1;
            P(:,j) = x.*(x.^2 - 3*y.^2);
            j = j + 1;
        case 4
            P(:,j) = x.*y.*(x.^2 - y.^2);
            j = j + 1;
            P(:,j) = y.*z.*(y.^2 - 3*x.^2);
            j = j + 1;
            P(:,j) = x.*y.*(-r2 + 7*z.^2);
            j = j + 1;
            P(:,j) = y.*z.*(3*r2 - 7*z.^2);
            j = j + 1;
            P(:,j) = 35*z.^4 - 30*r2.*z.^2 + 3*r2^2;
            j = j + 1;
            P(:,j) = x.*z.*(3*r2 - 7*z.^2);
            j = j + 1;
            P(:,j) = (x.^2 - y.^2).*(-r2 + 7*z.^2);
            j = j + 1;
            P(:,j) = x.*z.*(x.^2 - 3*y.^2);
            j = j + 1;
            P(:,j) = x.^4 - 6*x.^2.*y.^2 + y.^4;
            j = j + 1;
        case 5
            P(:,j) = y.*(5*x.^4-10*x.^2.*y.^2+ y.^4);
            j = j+1;
            P(:,j) = x.*y.*(x.^2-y.^2).*z;
            j = j+1;
            P(:,j) = y.*(-3*x.^2+y.^2).*(-r2+9*z.^2);
            j = j+1;
            P(:,j) = x.*y.*z.*(-r2+3*z.^2);
            j = j+1;
            P(:,j) = y.*(r2^2-14*r2*z.^2+21*z.^4);
            j = j+1;
            P(:,j) = z.*(15*r2^2-70*r2*z.^2+63*z.^4);
            j = j+1;
            P(:,j) = x.*(r2^2-14*r2*z.^2+21*z.^4);
            j = j+1;
            P(:,j) = (x.^2-y.^2).*z.*(-r2+3*z.^2);
            j = j+1;
            P(:,j) = x.*(x.^2-3*y.^2).*(-r2+9*z.^2);
            j = j+1;
            P(:,j) = (x.^4-6*x.^2.*y.^2+y.^4).*z;
            j = j+1;
            P(:,j) = x.*(x.^4-10*x.^2.*y.^2+5*y.^4);
            j = j+1;
        case 6
            P(:,j) = x.*y.*(3*x.^4-10*x.^2.*y.^2+3*y.^4);            
            j = j+1;
            P(:,j) = y.*(5*x.^4-10*x.^2.*y.^2+y.^4).*z;
            j = j+1;
            P(:,j) = x.*y.*(x.^2-y.^2).*(-r2+11*z.^2);
            j = j+1;
            P(:,j) = y.*(-3*x.^2+y.^2).*z.*(-3*r2+11*z.^2);
            j = j+1;
            P(:,j) = x.*y.*(r2-18*z.^2+33*z.^4);
            j = j+1;
            P(:,j) = y.*z.*(5*r2-30*z.^2+33*z.^4);
            j = j+1;
            P(:,j) = (-5*r2+105*z.^2-315*z.^4+231*z.^6);
            j = j+1;
            P(:,j) = x.*z.*(5*r2-30*z.^2+33*z.^4);
            j = j+1;
            P(:,j) = (x.^2-y.^2).*(1*r2-18*z.^2+33*z.^4);
            j = j+1;
            P(:,j) = x.*(x.^2-3*y.^2).*z.*(-3*r2+11*z.^2);
            j = j+1;
            P(:,j) = (x.^4-6*x.^2.*y.^2+y.^4).*(-r2+11*z.^2);
            j = j+1;
            P(:,j) = x.*(x.^4-10*x.^2.*y.^2+5*y.^4).*z;
            j = j+1;
            P(:,j) = (x.^6-15*x.^4.*y.^2+15*x.^2.*y.^4-y.^6);
            j = j+1;
        otherwise
            error('Only spherical harmonics up to degree 6 are currently supported.')
    end
end

end            
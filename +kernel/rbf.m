classdef (Abstract) rbf
%RBF Abstract class for radial basis function (RBF) kernels
%   Abstract class for implementing radial basis function (RBF) kernels for
%   various types.  In general a kernel has the form Phi(x,y), where x and y are
%   elements of R^d. A RBF kernel is one where Phi(x,y) = phi(r), where r=|x-y|
%   is the Euclidean distance between x and y. All RBF kernels should be
%   subclasses of this abstract class and impmlement functions for evaluating
%   the kernels at r, as well as the derivatives of the kernel with respect to
%   r. Specifically, two derivatives should be implemented for the rbf kernel
%   defined by phi(r):
%   eta(r) = (1/r) (d/dr phi(r))
%   zeta(r) = (1/r) (d/dr eta(r))
%   Together with phi(r), these functions are enough to define several second
%   order differential operators.
%
%   See also PHSEVEN, PHSODD, GAUSSIAN, IMQ, MQ

% Copyright 2024 by Grady B. Wright

    methods (Abstract)
        % Evaluation of the RBF phi(r)
        p = phi(rbfobj,r)
        % Evaluation of (1/r)(d/dr phi)
        d1p = eta(rbfobj,r)
        % Evaluation of (1/r)(d/dr eta)
        d2p = zeta(rbfobj,r)
    end

    methods
        function p = feval(obj,r)
            %PHI Evaluation of the kernel at r
            p = obj.phi(r);
        end
    end
end
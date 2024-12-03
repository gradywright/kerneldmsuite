classdef phsodd < kernel.rbf
%PHODD Polyharmonic spline radial basis function (RBF) kernel for odd dimensions
%   Implementation of the polyharmonic spline RBF kernel of order k in odd dimensions,
%   which (up to scaling) is defined as
%   phi(r) = (-1)^k*r^(2k-1)
%   Here k should be a positive integer.
%
%   To construct an odd polyharmonic spline RBF object use:
%   p = phsodd(k);
%
%   See also PHSEVEN and RBF

% Copyright 2024 by Grady B. Wright

    properties (Access=public)
        % Order of the polyharmonic spline kernel
        order
        sgn
    end

    methods (Access=public)

        function obj = phsodd(order)
            obj.order = order;
            obj.sgn = (-1)^(order);
        end

        function p = phi(obj,r)
            %PHI Evaluation of the polyharmonic spline kernel at r
            ell = obj.order;
            p = obj.sgn*(r.^(2*ell-1));
        end

        function d1p = eta(obj,r)
            %ETA Evaluation of (1/r)(d/dr phi(r))
            ell = obj.order;
            d1p = (obj.sgn*(2*ell-1))*(r.^(2*ell-3));
        end

        function d2p = zeta(obj,r)
            %ZETA Evaluation of (1/r)(d/dr eta(r))
            ell = obj.order;
            d2p = (obj.sgn*(3+4*ell*(ell-2)))*(r.^(2*ell-5));
            % When using this to compute an approximation to the Laplacian
            % the result will actually be defined at r=0 even for ell=2.
            % However, zeta will be undefined there.  The reason is that there
            % is and r^2 multiplying zeta in the Laplacian, so that the overall
            % result will be zero when r=0. We will fix zeta to be zero when r=0
            % to avoid these issues
            if ell == 2
                d2p(isinf(d2p)) = 0;
            end
        end

        % function out = subsref(p, index)
        %     idx = index(1).subs;
        %     switch index(1).type
        %         case '()'
        %             out = feval(p,idx{1});
        %         case '.'
        %             if strcmp(idx,'order')
        %                 out = p.order;
        %             elseif strcmp(idx,'sgn')
        %                 out = p.sgn;
        %             else
        %                 error('KDMSUITE:rbf:subsref','No such property')
        %             end
        %     end
        % end
    
    end
end
classdef phseven < kernel.rbf
%PHSEVEN Polyharmonic spline radial basis function (RBF) kernel for even dimensions
%   Implementation of the polyharmonic spline RBF kernel of order k in even dimensions,
%   which (up to scaling) is defined as
%   phi(r) = (-1)^k*r^(2k)*log(r)
%   Here k should be a positive integer.
%
%   To construct an even polyharmonic spline RBF object use:
%   p = phseven(k);
%
%   See also PHSODD

% Copyright 2024 by Grady B. Wright

    properties (Access=public)
        % Order of the polyharmonic spline kernel
        order
        sgn
    end

    methods (Access=public)

        function obj = phseven(order)
            obj.order = order;
            obj.sgn = (-1)^(order);
        end

        function p = phi(obj,r)
            %PHI Evaluation of the polyharmonic spline kernel at r
            ell = obj.order;
            % Handle the removable singularity at the origin
            r(r==0) = 1;
            p = obj.sgn*(r.^(2*ell).*log(r));
        end

        function d1p = eta(obj,r)
            %ETA Evaluation of (1/r)(d/dr phi(r))
            ell = obj.order;
            % Do in pieces to properly handle the removable singularity at the origin
            d1p = r.^(2*ell-2);
            % Deal with the removable singularity at the origin
            r(r==0) = 1;
            d1p = d1p + (2*ell)*(r.^(2*ell-2).*log(r));
            d1p = obj.sgn*d1p;

        end

        function d2p = zeta(obj,r)
            %ZETA Evaluation of (1/r)(d/dr eta(r))
            ell = obj.order;
            % Do in pieces to properly handle the removable singularity at the origin
            d2p = r.^(2*ell-4).*(-2 + 4*ell);
            % Deal with the removable singularity at the origin
            r(r==0) = 1;
            d2p = d2p + (4*ell*(ell-1))*(r.^(2*ell-4).*log(r));
            d2p = obj.sgn*d2p;
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
classdef imq < kernel.rbf
%IMQ Inverse multiquadric radial basis function (RBF) kernel
%   Implementation of the Inverse multiquadric RBF kernel, which is defined as
%   phi(r) = 1/sqrt(1 + (epsilon*r)^2)
%   where epsilon > 0 is the shape parameter.
%
%   To construct an inverse multiquadric RBF object use:
%   p = imq(epsilon);
%
%   See also RBF

% Copyright 2024 by Grady B. Wright

    properties (Access=public)
        % Shape parameter
        ep
    end

    properties (Access=private)
        % Square of the shape parameter
        epsq
        epsqsq
    end

    methods (Access=public)

        function obj = imq(ep)
            obj.ep = ep;
            obj.epsq = ep^2;
            obj.epsqsq = ep^4;
        end

        function p = phi(obj,r)
            %PHI Evaluation of the inverse multiquadric kernel at r
            ep2 = obj.epsq;
            p = 1./sqrt(1 + ep2*r.^2);
        end

        function d1p = eta(obj,r)
            %ETA Evaluation of (1/r)(d/dr phi(r))
            ep2 = obj.epsq;
            d1p = ep2./(1 + ep2*r.^2).^(3/2);
        end

        function d2p = zeta(obj,r)
            %ZETA Evaluation of (1/r)(d/dr eta(r))
            ep2 = obj.epsq;
            ep4 = obj.epsqsq;
            d2p = (-3*ep4)./(1 + ep2*r.^2).^(5/2);
        end

        % function out = subsref(p, index)
        %     idx = index(1).subs;
        %     switch index(1).type
        %         case '()'
        %             out = feval(p,idx{1});
        %         case '.'
        %             if strcmp(idx,'epsq')
        %                 out = p.epsq;
        %             elseif strcmp(idx,'epsq')
        %                 out = p.epsq;
        %             elseif strcmp(idx,'epsq')
        %                 out = p.epsq;
        %             else
        %                 error('KDMSUITE:rbf:subsref','No such property')
        %             end
        %     end
        % end
    end
end
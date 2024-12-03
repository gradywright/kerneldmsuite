classdef gaussian < kernel.rbf
%GAUSSIAN Gaussian radial basis function (RBF) kernel
%   Implementation of the Gaussian RBF kernel, which is defined as
%   phi(r) = exp(-(epsilon*r)^2
%   where epsilon > 0 is the shape parameter.
%
%   To construct a Gaussian RBF object use:
%   p = gaussian(epsilon);
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

        function obj = gaussian(ep)
            obj.ep = ep;
            obj.epsq = ep^2;
        end

        function p = phi(obj,r)
            %PHI Evaluation of the Gaussian kernel at r
            ep2 = obj.epsq;
            p = exp(-ep2*r.^2);
        end

        function d1p = eta(obj,r)
            %ETA Evaluation of (1/r)(d/dr phi(r))
            ep2 = obj.epsq;
            d1p = -2*ep2*exp(-ep2*r.^2);
        end

        function d2p = zeta(obj,r)
            %ZETA Evaluation of (1/r)(d/dr eta(r))
            ep2 = obj.epsq;
            d2p = (2*ep2)^2*exp(-ep2*r.^2);
        end

        % function out = subsref(p, index)
        %     idx = index(1).subs;
        %     switch index(1).type
        %         case '()'
        %             out = feval(p,idx{1});
        %         case '.'
        %             if strcmp(idx,'ep')
        %                 out = p.ep;
        %             elseif strcmp(idx,'epsq')
        %                 out = p.epsq;
        %             else
        %                 error('KDMSUITE:rbf:subsref','No such property')
        %             end
        %     end
        % end        
    end
end
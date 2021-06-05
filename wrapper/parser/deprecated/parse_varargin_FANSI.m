%% [mu1,alpha1,tol,maxiter,wmap,solver,constraint,b0dir]=parse_varargin_FANSI(arg)
%
% Description: parser for qsmFANSI.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 September 2017
% Date modified: 26 September 2017
% Date modified: 1 April 2019
% Date last modified: 27 Feb 2020 (v0.8.0)
%
% function [mu1,mu2,alpha1,tol,maxiter,wmap,solver,constraint,b0dir]=parse_varargin_FANSI(arg)
function [mu1,alpha1,wmap,options,b0dir]=parse_varargin_FANSI(arg)

% predefine parameters
alpha1  = 3e-5;
mu1     = 5e-5;
wmap    = [];
% maxiter = 40;
% solver = 'nonlinear';
% constraint = 'TGV';
% tol = 1;
b0dir = [0,0,1];
options.maxOuterIter    = 50;
options.tol_update      = 1;
options.tgv             = false;
options.nonlinear       = true;
options.isWeakHarmonic 	= false;
options.beta            = 1;
options.muh             = 1;
options.mu2             = 1;

% use user defined input if any
if ~isempty(arg)
    for kvar = 1:length(arg)
        if strcmpi(arg{kvar},'lambda')
            alpha1 = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'mu')
            mu1 = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'mu2')
            options.mu2 = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'tol')
%             tol = arg{kvar+1};
            options.tol_update = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'iteration')
%             maxiter = arg{kvar+1};
            options.maxOuterIter = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'weight')
            wmap = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'linear')
%             solver = 'linear';
            options.nonlinear = false;
        end
        if strcmpi(arg{kvar},'nonlinear')
%             solver = 'nonlinear';
            options.nonlinear = true;
        end
        if strcmpi(arg{kvar},'tv')
%             constraint = 'TV';
            options.tgv = false;
        end
        if strcmpi(arg{kvar},'tgv')
%             constraint = 'TGV';
            options.tgv = true;
        end
        if strcmpi(arg{kvar},'b0dir')
            b0dir = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'gradient_mode')
            gradient_mode = arg{kvar+1};
            switch gradient_mode
                case 'Vector field'
                    options.gradient_mode = 0;
                case 'L1 norm'
                    options.gradient_mode = 1;
                case 'L2 norm'
                    options.gradient_mode = 2;
            end
        end
        if strcmpi(arg{kvar},'isWeakHarmonic')
            options.isWeakHarmonic = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'beta')
            options.beta = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'muh')
            options.muh = arg{kvar+1};
        end
    end
end

end
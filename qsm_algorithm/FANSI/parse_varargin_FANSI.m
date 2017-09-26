%% [mu1,alpha1,tol,maxiter,wmap,solver,constraint,b0dir]=parse_varargin_FANSI(arg)
%
% Description: parser for qsmFANSI.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 September 2017
% Date last modified: 26 September 2017
%
function [mu1,alpha1,tol,maxiter,wmap,solver,constraint,b0dir]=parse_varargin_FANSI(arg)
alpha1 = 3e-5;
mu1 = 5e-5;
maxiter = 40;
wmap = [];
solver = 'nonlinear';
constraint = 'TGV';
tol = 1;
b0dir = [0,0,1];

if ~isempty(arg)
    for kvar = 1:length(arg)
        if strcmpi(arg{kvar},'lambda')
            alpha1 = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'mu')
            mu1 = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'tol')
            tol = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'iteration')
            maxiter = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'weight')
            wmap = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'linear')
            solver = 'linear';
        end
        if strcmpi(arg{kvar},'nonlinear')
            solver = 'nonlinear';
        end
        if strcmpi(arg{kvar},'tv')
            constraint = 'TV';
        end
        if strcmpi(arg{kvar},'tgv')
            constraint = 'TGV';
        end
        if strcmpi(arg{kvar},'b0dir')
            b0dir = arg{kvar+1};
        end
    end
end
end
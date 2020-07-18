%% function [lambda,tol,maxiter,wmap,initGuess,optimise,b0dir] = parse_varargin_iLSQR(arg)
%
% Description: parser for qsmIterativeLSQR.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 September 2017
% Date last modified: 26 September 2017
% Date last modified: 27 Feb 2020 (v0.8.0)
%
function [lambda,tol,maxiter,wmap,initGuess,optimise,b0dir] = parse_varargin_iLSQR(arg)

% predefine parameters
lambda      = 1e-1;
tol         = 1e-3;
maxiter     = 50;
wmap        = [];
initGuess   = [];
optimise    = false;
b0dir       = [0,0,1];

% use user defined input if any
if ~isempty(arg)
    for kvar = 1:length(arg)
        if strcmpi(arg{kvar},'lambda')
            lambda = arg{kvar+1};
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
        if strcmpi(arg{kvar},'initGuess')
            initGuess = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'optimise')
            optimise = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'b0dir')
            b0dir = arg{kvar+1};
        end
    end
end

end
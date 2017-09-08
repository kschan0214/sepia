%% function [lambda, tol, maxiter, wmap, initGuess, optimise] = parse_varargin_iLSQR(arg)
%
% Description: parser for qsmIterativeLSQR.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 September 2017
% Date last modified: 
%
function [lambda, tol, maxiter, wmap, initGuess, optimise] = parse_varargin_iLSQR(arg)
lambda = 1e-1;
tol = 1e-3;
maxiter = 50;
wmap = [];
initGuess = [];
optimise = false;

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
    end
end
end
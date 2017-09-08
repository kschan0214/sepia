%% function [lambda, optimise] = parse_varargin_CFL2norm(arg)
%
% Description: parser for qsmClosedFormL2.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 September 2017
% Date last modified: 
%
function [lambda, optimise] = parse_varargin_CFL2norm(arg)
lambda = 1e-1;
optimise = false;
if ~isempty(arg)
    for kvar = 1:length(arg)
        if strcmpi(arg{kvar},'lambda')
            lambda = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'optimise')
            optimise = arg{kvar+1};
        end
    end
end
end
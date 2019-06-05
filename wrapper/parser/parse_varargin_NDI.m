%% function [tolerance,stepSize,iteration,weight,b0dir] = parse_varargin_NDI(arg)
%
% Description: parser for NDI.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 5 June 2019
% Date last modified: 
%
function [tolerance,stepSize,iteration,weight,b0dir] = parse_varargin_NDI(arg)
% predefine parameters
tolerance = 1;
stepSize = 1;
iteration = 200;
b0dir=[0,0,1];
weight = [];
if ~isempty(arg)
    for kvar = 1:length(arg)
        if strcmpi(arg{kvar},'tol')
            tolerance = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'stepsize')
            stepSize = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'iteration')
            iteration = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'weight')
            weight = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'b0dir')
            b0dir = arg{kvar+1};
        end
    end
end
end
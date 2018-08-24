%% function [mask] = parse_varargin_UnwrapPhase_3DBestPath(arg)
%
% Description: parser for UnwrapPhaseMacro.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 8 September 2017
% Date last modified: 
%
function [mask] = parse_varargin_UnwrapPhase_3DBestPath(arg)
mask = [];
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'mask')
        mask = arg{kkvar+1};
        continue
    end
end
end
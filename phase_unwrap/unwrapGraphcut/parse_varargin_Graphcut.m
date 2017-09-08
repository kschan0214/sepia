%% [magn, subsampling] = parse_varargin_Graphcut(arg)
%
% Description: parser for UnwrapPhaseMacro.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 8 September 2017
% Date last modified: 
%
function [magn, subsampling] = parse_varargin_Graphcut(arg)
magn = [];
subsampling = 1;
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'Magn')
        magn = arg{kkvar+1};
        continue
    end
    if  strcmpi(arg{kkvar},'Subsampling')
        subsampling = arg{kkvar+1};
        continue
    end
end
end
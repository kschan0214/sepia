%% function [magn] = parse_varargin_RegionGrowing(arg)
%
% Description: parser for UnwrapPhaseMacro.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 8 September 2017
% Date last modified: 
%
function [magn] = parse_varargin_RegionGrowing(arg)
magn = [];
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'Magn')
        magn = arg{kkvar+1};
        continue
    end
end
end
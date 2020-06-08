%% [refine] = parse_varargin_VSHARPSTI(arg)
%
% Description: parser for BackgroundRemovalMacro.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 15 April 2017
% Date last modified: 
%
function [radius] = parse_varargin_VSHARPSTI(arg)
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'radius')
        radius = arg{kkvar+1};
        continue
    end
end
end
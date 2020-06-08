%% [radius, threshold, refine] = parse_varargin_SHARP(arg)
%
% Description: parser for BackgroundRemovalMacro.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 8 September 2017
% Date last modified: 
%
function [radius, threshold] = parse_varargin_SHARP(arg)
radius = 4;
threshold = 0.03;
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'radius')
        radius = arg{kkvar+1};
        continue
    end
    if  strcmpi(arg{kkvar},'threshold')
        threshold = arg{kkvar+1};
        continue
    end
end
end
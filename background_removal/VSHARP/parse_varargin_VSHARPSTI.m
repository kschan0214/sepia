%% [refine] = parse_varargin_VSHARPSTI(arg)
%
% Description: parser for BackgroundRemovalMacro.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 8 September 2017
% Date last modified: 
%
function [refine] = parse_varargin_VSHARPSTI(arg)
refine = false;
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'refine')
        refine = arg{kkvar+1};
        continue
    end
end
end
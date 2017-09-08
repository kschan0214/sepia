%% [radius, threshold, refine] = parse_varargin_RESHARP(arg)
%
% Description: parser for BackgroundRemovalMacro.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 8 September 2017
% Date last modified: 
%
function [radius, alpha, refine] = parse_varargin_RESHARP(arg)
radius = 4;
alpha = 0.01;
refine = false;
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'radius')
        radius = arg{kkvar+1};
        continue
    end
    if  strcmpi(arg{kkvar},'alpha')
        alpha = arg{kkvar+1};
        continue
    end
    if strcmpi(arg{kkvar},'refine')
        refine = arg{kkvar+1};
        continue
    end
end
end
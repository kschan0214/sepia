%% [iteration, refine] = parse_varargin_iHARPERELLA(arg)
%
% Description: parser for BackgroundRemovalMacro.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 8 September 2017
% Date last modified: 
%
function [iteration, refine] = parse_varargin_iHARPERELLA(arg)
iteration = 100;
refine = false;
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'iteration')
        iteration = arg{kkvar+1};
        continue
    end
    if strcmpi(arg{kkvar},'refine')
        refine = arg{kkvar+1};
        continue
    end
end
end
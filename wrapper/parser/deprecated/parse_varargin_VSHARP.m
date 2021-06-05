%% function [radius,refine] = parse_varargin_VSHARP(arg)
%
% Description: parser for BackgroundRemovalMacro.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 8 September 2017
% Date last modified: 
%
function [radius] = parse_varargin_VSHARP(arg)
radius = 5:-1:1;
if ~isempty(arg)
    for kvar = 1:length(arg)
        if strcmpi(arg{kvar},'radius')
            tmp = arg{kvar+1};
            radius = sort(tmp,'descend');
        end
    end
end
end
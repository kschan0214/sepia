%% function [thre_tkd,b0dir] = parse_varargin_Star(arg)
%
% Description: parser for qsmTKD.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 September 2017
% Date last modified: 27 Feb 2020 (v0.8.0)
%
function [padSize] = parse_varargin_Star(arg)

% predefine parameters
padSize = [12,12,12];

% use user defined input if any
if ~isempty(arg)
    for kvar = 1:length(arg)
        if strcmpi(arg{kvar},'padsize')
        	padSize = double(arg{kvar+1});
        end
    end
end

end
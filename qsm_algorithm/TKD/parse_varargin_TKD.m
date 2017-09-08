%% function thre_tkd = parse_varargin_TKD(arg)
%
% Description: parser for qsmTKD.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 September 2017
% Date last modified: 
%
function thre_tkd = parse_varargin_TKD(arg)
% predefine parameters
thre_tkd = 0.15;
if ~isempty(arg)
    for kvar = 1:length(arg)
        if strcmpi(arg{kvar},'threshold')
            thre_tkd = arg{kvar+1};
        end
    end
end
end
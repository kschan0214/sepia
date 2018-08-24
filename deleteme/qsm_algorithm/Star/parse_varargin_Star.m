%% function [thre_tkd,b0dir] = parse_varargin_Star(arg)
%
% Description: parser for qsmTKD.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 September 2017
% Date last modified: 
%
function [te,padSize,b0,b0dir] = parse_varargin_Star(arg)
% predefine parameters
te = 40;
b0dir=[0,0,1];

if ~isempty(arg)
    for kvar = 1:length(arg)
        if strcmpi(arg{kvar},'te')
            te = arg{kvar+1}*1e3;
        end
        if strcmpi(arg{kvar},'b0dir')
            b0dir = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'TE')
            te = double(arg{kvar+1}) * 1e3;
        end
        if strcmpi(arg{kvar},'B0')
            b0 = double(arg{kvar+1});
        end
        if strcmpi(arg{kvar},'padsize')
        	padSize = double(arg{kvar+1});
        end
    end
end
end
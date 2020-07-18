%% function [b0,b0dir,te] = parse_varargin_QSMmacro(arg)
%
% Description: parser for QSMMacro.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 27 Feb 2020
% Date last modified: 
%
function [b0,b0dir,te,mask_ref] = parse_varargin_QSMmacro(arg)

% predefine parameters
b0                  = 3; % T
b0dir               = [0,0,1];
te                  = 40e-3; % s
mask_ref            = []; 

if ~isempty(arg)
    for kvar = 1:length(arg)
        if strcmpi(arg{kvar},'B0')
            b0 = double(arg{kvar+1});
        end
        if strcmpi(arg{kvar},'b0dir')
            b0dir = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'TE')
            te = double(arg{kvar+1});
        end
        if strcmpi(arg{kvar},'reference_mask')
            mask_ref = arg{kvar+1};
        end
    end
end

end
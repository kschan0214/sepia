%% [B0_dir, tol, iteration, CGdefault, N_std, refine] = parse_varargin_PDF(matrixSize,arg)
%
% Description: parser for BackgroundRemovalMacro.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 September 2017
% Date last modified: 
%
function [tol, depth, peel] = parse_varargin_LBV(arg)
tol = 1e-4;
depth = 4;
peel = 1;
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'tol')
        tol = arg{kkvar+1};
        continue
    end
    if  strcmpi(arg{kkvar},'depth')
        depth = arg{kkvar+1};
        continue
    end
    if strcmpi(arg{kkvar},'peel')
        peel = arg{kkvar+1};
        continue
    end
end
end
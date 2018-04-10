%% [B0_dir, tol, iteration, CGdefault, N_std, refine] = parse_varargin_PDF(matrixSize,arg)
%
% Description: parser for qsmFANSI.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 September 2017
% Date last modified: 
%
function [B0_dir, tol, iteration, padSize, N_std, refine] = parse_varargin_PDF(arg)
B0_dir = [0,0,1];
tol = 0.1;
iteration = 30;
padSize = 40;
N_std = [];
refine = false;
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'b0dir')
        B0_dir = arg{kkvar+1};
        continue
    end
    if  strcmpi(arg{kkvar},'tol')
        tol = arg{kkvar+1};
        continue
    end
    if strcmpi(arg{kkvar},'iteration')
        iteration = arg{kkvar+1};
        continue
    end
    if strcmpi(arg{kkvar},'pad')
        padSize = arg{kkvar+1};
        continue
    end
    if strcmpi(arg{kkvar},'noisestd')
        N_std = arg{kkvar+1};
        continue
    end
    if strcmpi(arg{kkvar},'refine')
        refine = arg{kkvar+1};
        continue
    end
end
end

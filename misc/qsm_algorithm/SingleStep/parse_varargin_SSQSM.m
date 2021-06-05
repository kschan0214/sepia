%% [lambda,magn,tol,maxiter,Kernel_Sizes]=parse_varargin_SSQSM(arg)
%
% Description: parser for qsmFANSI.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 September 2017
% Date last modified: 26 September 2017
%
function [lambda,magn,tol,maxiter,Kernel_Sizes,b0dir]=parse_varargin_SSQSM(arg)
% function [B0,TE,lambda,magn,tol,maxiter,Kernel_Sizes]=parse_vararginSSQSM(arg)
% B0 = 3;
% TE = 1;             %second
lambda = 2.9e-2;
magn = [];
maxiter = 30;
tol = 1e-2;
Kernel_Sizes = 11:-2:3;
b0dir=[0,0,1];

if ~isempty(arg)
    for kvar = 1:length(arg)
%         if strcmpi(arg{kvar},'fieldStrength')
%             B0 = arg{kvar+1};
%         end
%         if strcmpi(arg{kvar},'te')
%             TE = arg{kvar+1};
%         end
        if strcmpi(arg{kvar},'tol')
            tol = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'iteration')
            maxiter = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'magnitude')
            magn = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'lambda')
            lambda = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'vkernel')
            if ~isempty(arg{kvar+1})
                Kernel_Sizes = arg{kvar+1};
            end
        end
        if strcmpi(arg{kvar},'b0dir')
            b0dir = arg{kvar+1};
        end
    end
end
end
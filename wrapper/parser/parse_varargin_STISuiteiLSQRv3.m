%% params = parse_varargin_STISuiteiLSQR(arg)
%
% Description: parser for QSM_iLSQR.p
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 September 2017
% Date last modified: 
%
% function [te,b0,b0dir,padSize,VoxelSize] = parse_varargin_STISuiteiLSQRv3(arg)
function params = parse_varargin_STISuiteiLSQRv3(arg)
% params.H=[0 0 1];
params.niter=100;
% params.TE=40;
% params.B0=3;
params.tol_step1=0.01;
params.tol_step2=0.001;
params.Kthreshold=0.25;
params.padsize=[4,4,4];

for kvar = 1:length(arg)
%     if strcmpi(arg{kvar},'b0dir')
%         params.H = double(arg{kvar+1});
%         continue
%     end
    if strcmpi(arg{kvar},'tol_step1')
        params.tol_step1 = double(arg{kvar+1});
        continue
    end
    if strcmpi(arg{kvar},'tol_step2')
        params.tol_step2 = double(arg{kvar+1});
        continue
    end
    if strcmpi(arg{kvar},'iteration')
        params.niter = double(arg{kvar+1});
        continue
    end
%     if strcmpi(arg{kvar},'TE')
%         params.TE = double(arg{kvar+1}) * 1e3;
%         continue
%     end
%     if strcmpi(arg{kvar},'B0')
%         params.B0 = double(arg{kvar+1});
%         continue
%     end
    if strcmpi(arg{kvar},'threshold')
        params.Kthreshold = double(arg{kvar+1});
        continue
    end
    if strcmpi(arg{kvar},'padsize')
        params.padsize = double(arg{kvar+1});
        continue
    end
end
end
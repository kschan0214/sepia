%% [b0,b0dir,voxelSize,matrixSize,te,dte,CF]=SyntheticQSMHubHeader(multiecho,varargin)
%
% Usage:
%
% Input
% --------------
%
% Output
% --------------
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 11 April 2018
% Date last modified:
%
%
function [b0,b0dir,voxelSize,matrixSize,te,dte,CF]=SyntheticQSMHubHeader(nii,varargin)
gyro = 42.57747892;

a=qGetR([0, nii.hdr.hist.quatern_b,nii.hdr.hist.quatern_c,nii.hdr.hist.quatern_d]);
B0_dir = -a(3,:);
voxelSize = nii.hdr.dime.pixdim(2:4);
matrixSize = nii.hdr.dime.dim(2:4);
nt = nii.hdr.dime.dim(5);

% [matrixSize(1),matrixSize(2),matrixSize(3),nt] = size(multiecho);

[b0,b0dir,voxelSize,matrixSize,te] = process_options(varargin,'b0',3,'b0dir',B0_dir,...
    'voxel',voxelSize,'matrixsize',matrixSize,'te',linspace(1e-3,30e-3,nt));

if length(te)<2
    dte = te;
else
    dte = te(2)-te(1); %s
end

CF = gyro * 1e6 * b0;
            
end
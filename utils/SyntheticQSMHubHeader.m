%% [b0,b0dir,voxelSize,matrixSize,te,dte,CF]=SyntheticQSMHubHeader(nii,varargin)
%
% Input
% --------------
% nii           : NIfTI matlab structure
% varargin (for user input)
% ---------
% 'b0'          : field strength in tesla (default: 3)
% 'b0dir'       : main magnetic field direction (default: based on quantern of NIfTI header)
% 'voxel'       : voxel size in mm
% 'matrixsize'  : size of the input image
% 'te'          : echo times in second
%
% Output
% --------------
% b0            : field strength in tesla (default: 3)
% b0dir         : main magnetic field direction
% voxelSize     : voxel size in mm
% matrixSize    : size of the input image
% te            : echo times in second
% dte           : echo spacing in second
% CF            : imaging frequency
%
% Description: This function generates a synthetic data header for
% qsm_hub.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 11 April 2018
% Date last modified: 19 April 2018
%
%
function [b0,b0dir,voxelSize,matrixSize,te,dte,CF]=SyntheticQSMHubHeader(nii,varargin)

gyro = 42.57747892;

% get main field direction
a=qGetR([0, nii.hdr.hist.quatern_b,nii.hdr.hist.quatern_c,nii.hdr.hist.quatern_d]);
B0_dir = -a(3,:);
% get voxel size
voxelSize = nii.hdr.dime.pixdim(2:4);
% get matrix size
matrixSize = nii.hdr.dime.dim(2:4);
% get no. of echoes
nt = nii.hdr.dime.dim(5);

% generate the variables needed 
[b0,b0dir,voxelSize,matrixSize,te] = process_options(varargin,'b0',3,'b0dir',B0_dir,...
    'voxel',voxelSize,'matrixsize',matrixSize,'te',linspace(1e-3,30e-3,nt));

% compute the echo spacing
if length(te)<2
    dte = te;
else
    dte = te(2)-te(1); %s
end

% compute the imaging frequency
CF = gyro * 1e6 * b0;
            
end
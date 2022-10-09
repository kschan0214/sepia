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
% Date modified: 24 May 2018
% Date modified: 21 Jan 2021
%
%
function [B0,B0_dir,voxelSize,matrixSize,TE,delta_TE,CF]=SyntheticQSMHubHeader(nii,varargin)

gyro = 42.57747892;

[B0,B0_dir,voxelSize,matrixSize,TE]=parse_varargin(varargin);

if isempty(B0_dir)
    % get main field direction
    B0_dir = get_B0_dir_from_nifti(nii);
%     a = sqrt(1 - nii.hdr.hist.quatern_b^2 - nii.hdr.hist.quatern_c^2 - nii.hdr.hist.quatern_d^2);
%     rotmat_nifti=qGetR([a, nii.hdr.hist.quatern_b,nii.hdr.hist.quatern_c,nii.hdr.hist.quatern_d]);
%     B0_dir = rotmat_nifti \ [0;0;1];
end
if isempty(voxelSize)
    % get voxel size
    voxelSize = nii.hdr.dime.pixdim(2:4);
end
if isempty(matrixSize)
    % get matrix size
    matrixSize = nii.hdr.dime.dim(2:4);
end
if isempty(TE)
    % get no. of echoes
    nt = nii.hdr.dime.dim(5);
    TE = linspace(1e-3,30e-3,nt);
end

% generate the variables needed 
% [b0,b0dir,voxelSize,matrixSize,te] = process_options(varargin,'b0',3,'b0dir',B0_dir,...
%     'voxel',voxelSize,'matrixsize',matrixSize,'te',linspace(1e-3,30e-3,nt));

% compute the echo spacing
if length(TE)<2
    delta_TE = TE;
else
    delta_TE = TE(2)-TE(1); %s
end

% compute the imaging frequency
CF = gyro * 1e6 * B0;
            
end

function [b0,b0dir,voxelSize,matrixSize,te]=parse_varargin(arg)
% set default header
b0 = 3; %T
b0dir = [];
voxelSize = [];
matrixSize = [];
te = [];

if ~isempty(arg)
    for kvar = 1:length(arg)
        if strcmpi(arg{kvar},'b0')
            b0 = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'b0dir')
            b0dir = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'voxel')
            voxelSize = arg{kvar+1};
        end
        if  strcmpi(arg{kvar},'matrixsize')
            matrixSize = arg{kvar+1};
        end
        if  strcmpi(arg{kvar},'te')
            te = arg{kvar+1};
        end
    end
end
end
% BET.m
%   [Mask] = BET(iMag,matrix_size,voxel_size)
% 
%   output
%   Mask - the brain mask
% 
%   input
%   iMag - magnitude image
%   matrix_size - dimension of the 3D image stack
%   voxel_size - dimensions of the voxels 
%

function Mask = BET(iMag,matrix_size,voxel_size)
    % data type conversion
    matrix_size = double(matrix_size);
    voxel_size = double(voxel_size);
    n_vox = matrix_size(1)*matrix_size(2)*matrix_size(3);
    fM = double(reshape(iMag, [1,n_vox]));

try
    % c++ interface
    tmp = bet2(fM,matrix_size,voxel_size);
    Mask = reshape(tmp,matrix_size);
catch
    [STR, NAM, EXT] = fileparts(mfilename('fullpath'));
    error(['bet2.' mexext ' not found. Please run ''runmex'' in ' STR])
end

end


% Laplacian Boundary Value (LBV)
%   [RDF] = LBV(iFreq,Mask,matrix_size,voxel_size,tol,depth,peel)
% 
%   output
%   RDF - the relative difference field, or local field
% 
%   input
%   iFreq - the unwrapped field map
%   N_std - the noise standard deviation on the field map. (1 over SNR for single echo)
%   Mask - a binary 3D matrix denoting the Region Of Interest
%   matrix_size - the size of the 3D matrix
%   voxel_size - the size of the voxel in mm
%   tol(optional) - tolerance level
%   depth(optional) - multigrid level
%   peel(optional) - thickness of the boundary layer to be peeled off
%
%   written by Dong Zhou
%   zhou.dong@gmail.com
%   created:   : 6.12.2013
%   last modify: 6.24.2013


function [fL] = LBV(iFreq,Mask,matrix_size,voxel_size,tol,depth,peel)
    if (nargin <6)
        depth = -1;
    end
    if (nargin < 7)
        peel = 0;
    end

    % data type conversion
    n_vox = matrix_size(1)*matrix_size(2)*matrix_size(3);
    fT = double(reshape(iFreq, [1,n_vox]));
    mask = double(reshape(Mask, [1,n_vox]));
    matrix_size = double(matrix_size);
    voxel_size = double(voxel_size);
    
    tmp = mexMGv6(fT,mask,matrix_size,voxel_size,tol,depth,peel,30,100,100);
    fL = reshape(tmp,matrix_size);

end

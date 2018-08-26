% LBV.m
%   [fL] = LBV(iFreq,Mask,matrix_size,voxel_size,tol,depth,peel,N1,N2,N3)
% 
%   output
%   fL - the local field
% 
%   input
%   iFreq - total magnetic field
%   Mask - ROI
%   matrix_size - dimension of the 3D image stack
%   voxel_size - dimensions of the voxels 
%   tol - iteration stopping criteria on the coarsest grid
%   depth - number of length scales. The largest length scale is 2^depth * voxel size.
%   peel - number of boundary layers to be peeled off
%   N1 - iterations on each depth before the recursive call
%   N2 - iterations on each depth after the recursive call
%   N3 - iterations on the finest scale after the FMG is finished.
%
%   When using the code, please cite 
%   Zhou et al. NMR in Biomed 27 (3), 312-319, 2014
% 
%   Created by Dong Zhou (zhou.dong@gmail.com) on 2013.06.12
%   Last modified by Dong Zhou on 2013.06.24

function [fL] = LBV(iFreq,Mask,matrix_size,voxel_size,tol,depth,peel,N1,N2,N3)
    if (nargin <5)
        tol = 0.01;
    end
    if (nargin <6)
        depth = -1;
    end
    if (nargin < 7)
        peel = 0;
    end
    if (nargin < 8)
        N1 = 30;
    end
    if (nargin < 9)
        N2 = 100;
    end
    if (nargin < 10)
        N3 = 100;
    end

    % data type conversion
    matrix_size = double(matrix_size);
    voxel_size = double(voxel_size);
    n_vox = matrix_size(1)*matrix_size(2)*matrix_size(3);
    fT = double(reshape(iFreq, [1,n_vox]));
    mask = double(reshape(Mask, [1,n_vox]));

    % c++ interface
    tmp = mexMGv6(fT,mask,matrix_size,voxel_size,tol,depth,peel,N1,N2,N3);
    fL = reshape(tmp,matrix_size);

end


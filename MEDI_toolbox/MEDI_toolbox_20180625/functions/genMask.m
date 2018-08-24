% Mask generation
%   Mask = genMask(iField, voxel_size)
% 
%   output
%   Mask - the biggest contiguous region that has decent SNR
%
%   input
%   iField - the complex MR image
%   voxel_size - the size of a voxel
%
%   Created by Tian Liu in 20013.07.24
%   Last modified by Tian Liu on 2013.07.24

function Mask = genMask(iField, voxel_size)
    matrix_size = size(iField(:,:,:,1));
    iMag = sqrt(sum(abs(iField).^2,4));
    m = iMag>(0.05*max(iMag(:)));           % simple threshold
    m1 = SMV(m,matrix_size, voxel_size, 10)>0.999;   % erode the boundary by 10mm
    l = bwlabeln(m1,6);                       % find the biggest contiguous region
    Mask = (l==mode(l(l~=0)));
%     Mask1 = SMV(Mask, matrix_size, voxel_size, 10)>0.001; % restore the enrosion by 4mm
end

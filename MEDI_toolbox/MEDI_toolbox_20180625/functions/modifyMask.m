% Modify the mask to exclude unreliable field points at the boundary
%   eliminate voxels with large field variations with 
%   their neighbors
%   these voxels are likely caused by unwrapping problems
%   or low signal intensity
%   [m] = modifyMask(mask,iFreq,vs,n)
%
%   output
%   m - new mask
% 
%   input
%   mask - ROI
%   iFreq - the field map
%   vs - voxel size
%   threshold for outliers is set as n * sigma
%
%   Created by Dong Zhou (zhou.dong@gmail.com) on 2013.07.20
%   Last modified by Dong Zhou on 2013.07.20


function [m] = modifyMask(mask,iFreq,vs,n)

if (nargin < 4) 
    n = 2;  % by default, 2sigma ~ 95%
end

[gx gy gz] = gradient(iFreq, vs(1),vs(2),vs(3));

sx = std(gx(:));
sy = std(gy(:));
sz = std(gz(:));

avg_gx = mean(gx(:))
avg_gy = mean(gy(:))
avg_gz = mean(gz(:))

m = (abs(gx-avg_gx)<sx*n) .* (abs(gy-avg_gy)<sy*n) .* (abs(gz-avg_gz)<sz*n);

m = m .* mask;

% pick the largest connected component
CC = bwconncomp(m,6);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);
m = zeros(size(m));
m(CC.PixelIdxList{idx}) = 1;

end






%% mask_r2s_morph = refine_brain_mask_using_r2s(r2s,mask,voxelSize,min_open,min_close)
%
% Input
% --------------
% r2s           : R2* map
% mask          : brain mask
% voxelSize     : spatial resolution of the input data
% min_open      : minium radius for imopen operation
% max_close     : minium radius for imclose operation
%
% Output
% --------------
% mask_r2s_morph: refined brain mask
%
% Description: refining brain mask by thresholding high R2* values on brain
% edge
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 4 October 2022
% Date last modified:
%
%
function mask_r2s_morph = refine_brain_mask_using_r2s(r2s,mask,voxelSize,min_open,min_close)

if nargin<5
    min_close = 2; % mm     % minimum imclode radius
end
if nargin<4
    min_open = 1; % mm      % minimum imopen radius
end

% excluding minimum from statistic
min_r2s = min(r2s(mask>0));

% compute stats
iqr_r2s     = iqr(r2s(and( mask>0, r2s>min_r2s )));
median_r2s  = median(r2s(and( mask>0, r2s>min_r2s )));

% assume everthing outside 3*IQR from median to be outliers
mask_r2s = r2s <= (median_r2s + 3*iqr_r2s) .* mask;

% basic morphology operation
open_radius_mm  = min(min_open, max(voxelSize)); open_radius_mm( open_radius_mm<1 ) = 1;
close_radius_mm = max(min_close, max(voxelSize));

% convert mm to voxel
open_radius_voxel   = min(round(open_radius_mm ./ voxelSize));
close_radius_voxel  = min(round(close_radius_mm ./ voxelSize));

% remove disconnected voxels
mask_r2s_morph = imopen(mask_r2s, strel('sphere', open_radius_voxel));
% get the largest single object
mask_r2s_morph = getLargestObject(mask_r2s_morph);
% reconnect voxels
mask_r2s_morph = imclose(mask_r2s_morph, strel('sphere', close_radius_voxel));
% make sure no holes in the centre of the brain
mask_r2s_morph = imfill(mask_r2s_morph,'holes');

end
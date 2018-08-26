function Mask_ROI_CSF = extract_CSF(R2s, Mask, voxel_size, flag_erode, thresh_R2s)

if isempty(R2s)
    Mask_ROI_CSF=[];
    return
end

if nargin < 5
    thresh_R2s = 5;
end
if nargin < 4
    flag_erode = 1;
end

n_region_cen = 3;

matrix_size = size(Mask);



% Center region (sphere)
[X,Y,Z] = ndgrid((1:matrix_size(1))*voxel_size(1),(1:matrix_size(2))*voxel_size(2),(1:matrix_size(3))*voxel_size(3));
X_cen = sum(X(Mask>0))/sum(Mask(:));
Y_cen = sum(Y(Mask>0))/sum(Mask(:));
Z_cen = sum(Z(Mask>0))/sum(Mask(:));

radius_cen = 30;
Mask_cen = sqrt(abs(X - X_cen).^2 + abs(Y - Y_cen).^2 + abs(Z - Z_cen).^2) <= radius_cen;


if flag_erode
    Mask = SMV(Mask, matrix_size, voxel_size, 10)>0.999;
end



Mask_raw_1 = (R2s < thresh_R2s).*Mask_cen;
CC = bwconncomp(Mask_raw_1,6);
numPixels = cellfun(@numel,CC.PixelIdxList);
[numPixels_sorted,idxs] = sort(numPixels,2,'descend');
ROIs_region_cen = zeros(matrix_size);
for i = 1:n_region_cen
    numPixels_sorted(i);
    idx = idxs(i);
    ROIs_region_cen(CC.PixelIdxList{idx}) = i;
end

Mask_raw_2 = (R2s < thresh_R2s).*Mask;
CC = bwconncomp(Mask_raw_2,6);
numPixels = cellfun(@numel,CC.PixelIdxList);
[numPixels_sorted,idxs] = sort(numPixels,2,'descend');
ROIs_region = zeros(matrix_size);
for i = 1:length(idxs)
    numPixels_sorted(i);
    idx = idxs(i);
    ROIs_region(CC.PixelIdxList{idx}) = i;
end

% Choose regions which appear at center
Mask_ROI_CSF = zeros(matrix_size);
for i = transpose(unique(ROIs_region(ROIs_region_cen>0 & ROIs_region > 0)))
    Mask_ROI_CSF(ROIs_region == i) = 1;
end

Mask_ROI_CSF(Mask == 0) = 0;

end
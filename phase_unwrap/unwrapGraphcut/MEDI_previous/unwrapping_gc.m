% unwrapping using graph cut with magnitude weighting 
% OUTPUT: unwph
%% How to use: 
% Unwrapping with full matrix: [unwph] = unwrapping_gc(iFreq_raw,iMag,voxel_size);
% Unwrapping with subsampling for speed up: [unwph] = unwrapping_gc(iFreq_raw,iMag,voxel_size,2);
% Jianwu Dong 2014/9/3

function [unwph, kk] = unwrapping_gc(iFreq_raw,iMag,voxel_size,SUBSAMPLE)


if nargin<4
    SUBSAMPLE = 1;
end

p = 2;
iFreq_raw(isnan(iFreq_raw))=0;
iFreq_raw1 = iFreq_raw;

if SUBSAMPLE == 2
    % subsample to speed up
    [sx sy sz] = size(iFreq_raw1);
    START = 1;
    allX = 1:sx;
    allY = 1:sy;
    allZ = 1:sz;
    subX = START:SUBSAMPLE:sx;
    subY = START:SUBSAMPLE:sy;
    subZ = START:SUBSAMPLE:sz;
    iFreq_raw2 = iFreq_raw1(START:SUBSAMPLE:end,START:SUBSAMPLE:end,:);
    iMag2 = iMag(START:SUBSAMPLE:end,START:SUBSAMPLE:end,:);
    voxel_size2(1) = voxel_size(1)/0.5;
    voxel_size2(2) = voxel_size(2)/0.5;
    voxel_size2(3) = voxel_size(3);
    [unwph2,iter,erglist] = phase_unwrap_3d(iFreq_raw2,p,iMag2,voxel_size2); 

   
    [ALLX,ALLY] = meshgrid(allY(:),allX(:));
    [SUBX,SUBY] = meshgrid(subY(:),subX(:));
    fm = zeros(size(iFreq_raw1));
    for ind = 1:sz
        fm(:,:,ind) = interp2(SUBX,SUBY,unwph2(:,:,ind),ALLX,ALLY,'*spline');
    end
    
    kk = round((fm - iFreq_raw1)/2/pi);
    unwph = iFreq_raw1 + 2*pi*kk;
    
else
    [unwph,iter,erglist] = phase_unwrap_3d(iFreq_raw1,p,iMag,voxel_size); 
end

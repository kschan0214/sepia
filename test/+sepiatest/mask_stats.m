%% stats = sepiatest.mask_stats(img, maskLogical)
%
% Description: computes a small set of summary statistics of img inside
% maskLogical, used as the regression-comparison payload instead of a
% full-resolution array diff (which is brittle across MATLAB/BLAS/FFT
% library versions for iterative solvers). Also returns a coarse
% block-mean downsample of the full array so gross spatial regressions
% (wrong sign, collapsed pattern, wrong scaling) are still caught.
%
% Input
% --------------
% img         : numeric array (any size)
% maskLogical : logical array, same size as img
%
% Output
% --------------
% stats : struct with fields mean, std, median, p1, p5, p25, p75, p95, p99,
%         nVoxelsInMask, nNaN, nInf (NaN/Inf counts inside the mask), and
%         downsampled (a small fixed-size block-mean array, independent of
%         img's native resolution)
%
function stats = mask_stats(img, maskLogical)

v = double(img(maskLogical));

stats = struct();
stats.nVoxelsInMask = numel(v);
stats.nNaN = sum(isnan(v));
stats.nInf = sum(isinf(v));

% percentile/summary stats computed on finite values only, so a handful of
% NaN/Inf voxels (already separately reported above) don't corrupt them
vFinite = v(isfinite(v));
if isempty(vFinite)
    stats.mean = NaN; stats.std = NaN; stats.median = NaN;
    stats.p1 = NaN; stats.p5 = NaN; stats.p25 = NaN;
    stats.p75 = NaN; stats.p95 = NaN; stats.p99 = NaN;
else
    stats.mean   = mean(vFinite);
    stats.std    = std(vFinite);
    stats.median = median(vFinite);
    pct = prctile(vFinite, [1 5 25 75 95 99]);
    stats.p1  = pct(1); stats.p5  = pct(2); stats.p25 = pct(3);
    stats.p75 = pct(4); stats.p95 = pct(5); stats.p99 = pct(6);
end

stats.downsampled = coarse_downsample(double(img) .* maskLogical, [8 8 6]);

end

%% block-mean downsample to a fixed target size (independent of native resolution)
function out = coarse_downsample(img, targetSize)

sz = size(img);
sz(end+1:3) = 1; % pad in case img is 2D
targetSize = min(targetSize, sz(1:3)); % never upsample

out = zeros(targetSize);
edges = cell(1,3);
for d = 1:3
    edges{d} = round(linspace(0, sz(d), targetSize(d)+1));
end

for ix = 1:targetSize(1)
    xr = (edges{1}(ix)+1):edges{1}(ix+1);
    for iy = 1:targetSize(2)
        yr = (edges{2}(iy)+1):edges{2}(iy+1);
        for iz = 1:targetSize(3)
            zr = (edges{3}(iz)+1):edges{3}(iz+1);
            block = img(xr,yr,zr);
            out(ix,iy,iz) = mean(block(:));
        end
    end
end

end

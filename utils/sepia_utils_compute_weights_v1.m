%% weights = sepia_utils_compute_weights_v1(fieldmapSD,mask)
%
% Input
% --------------
% fieldmapSD    : standard deviation of fieldmap
% mask          : signal mask
%
% Output
% --------------
% weights       : weighting map
%
% Description: Weights for QSM inversion, median
%
% Kwok-shing Chan @ DCCN
% kwokshing.chan@donders.ru.nl
% Date created: 21 Feb 2022
% Date modified:
%
%
function weights = sepia_utils_compute_weights_v1(fieldmapSD,mask)

if nargin < 2
    mask = ones(size(fieldmapSD));
end

weights                 = 1./fieldmapSD;
weights(isinf(weights)) = 0;
weights(isnan(weights)) = 0;

weights_1d = weights(mask>0);

% Step 1: robust z-score normalistion, centre at 0, S.D. of 1
iqr_mask    = iqr(weights_1d);
median_mask = median(weights_1d);
weights_1d = (weights_1d - median_mask) / iqr_mask;

% Step 2: re-scale and re-centre at 1  
weights_1d = (weights_1d/3 + 1);

% Clipping extreme values, assume anything outside 3*IQR is an outlier
% upper bound uses 3*IQR,  lower bound uses 0
iqr_mask    = iqr(weights_1d);      % update IQR
median_mask = median(weights_1d);   % update median

ub = median_mask + 3*iqr_mask;
lb = 0;
weights_1d(weights_1d > ub) = ub;
weights_1d(weights_1d < lb) = lb;

% normalise to [0,1]
weights_1d = (weights_1d - lb) / (ub - lb);

% put the normalised weights back to image
weights(mask>0) = weights_1d;

end
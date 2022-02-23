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

%% Step 1: Inversion of SD
weights                 = 1./fieldmapSD;
weights(isinf(weights)) = 0;
weights(isnan(weights)) = 0;

%  get all masked data
weights_1d = weights(mask>0);

%% Step 2: Normalisation
% get stats
iqr_mask    = iqr(weights_1d);
median_mask = median(weights_1d);

% the histogram of the weights is usually negatively skewed, so median +3*IQR should be safe
% for the right side of the distribution
weights_1d = weights_1d./(median_mask + 3*iqr_mask);
% weights_1d(weights_1d>1) = 1;
% % weights                 = weights./max(weights(mask>0));

%% Step 3: Re-centre
weights_1d = weights_1d - median(weights_1d) + 1;

% put the normalised weights back to image
weights(mask>0) = weights_1d;

%% TODO 20220223: could provide this as option in the future
% %% Step 4: Handle outliers on the right hand side of the histogram
% iqr_mask    = iqr(weights_1d);
% median_mask = median(weights_1d);
% threshold   = median_mask + 3*iqr_mask;
% 
% % use median to filter outliers
% % weights_filter = medfilt3(weights.*mask,[3, 3, 3]);
% weights_filter = imgaussfilt3(weights.*mask);
% % avoid introducing zero on the boundary, better to be a bit conservative here
% weights_filter(weights_filter == 0) = weights(weights_filter == 0);
% weights_final = weights;
% % replace outliers by median
% weights_final(weights>threshold) = weights_filter(weights>threshold);
% 
% weights = weights_final;

end
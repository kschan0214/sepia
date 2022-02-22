%% weights = sepia_utils_compute_weights_v0p8(fieldmapSD,mask)
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
% Description: (For legacy purpose) This is how SEPIA calculate the weights
% before v1.0
%
% Kwok-shing Chan @ DCCN
% kwokshing.chan@donders.ru.nl
% Date created: 21 Feb 2022
% Date modified:
%
%
function weights = sepia_utils_compute_weights_v0p8(fieldmapSD,mask)

if nargin < 2
    mask = ones(size(fieldmapSD));
end

weights                 = 1./fieldmapSD;
weights(isinf(weights)) = 0;
weights(isnan(weights)) = 0;
weights                 = weights./max(weights(mask>0));
    
end
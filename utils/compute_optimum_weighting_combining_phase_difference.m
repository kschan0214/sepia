%% [weights, fieldmapSD]= compute_optimum_weighting_combining_phase_difference(magn, TE)
%
% Input
% --------------
% magn          : multi-echo magnitude images, 4D
% TE            : echo times
%
% Output
% --------------
% weights       : weighting factors, 4D
% fieldmapSD    : theorectical noise derived from magnitude images
%
% Description: Robinson et al. 2017 NMR Biomed Appendix A2
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 13 August 2021 (v1.0)
% Date modified:
%
%
function [weights, fieldmapSD]= compute_optimum_weighting_combining_phase_difference(magn, TE)

dims    = size(magn);
dims(4) = dims(4) - 1;
    
fieldmapSD = zeros(dims, 'like', magn);   % fieldmap SD
for k = 1:dims(4)
    fieldmapSD(:,:,:,k) = 1./(TE(k+1)-TE(1)) ...
        * sqrt((magn(:,:,:,1).^2+magn(:,:,:,k+1).^2)./((magn(:,:,:,1).*magn(:,:,:,k+1)).^2));
end

% weights are inverse of the field map variance
weights = bsxfun(@rdivide,1./(fieldmapSD.^2),sum(1./(fieldmapSD.^2),4));

end
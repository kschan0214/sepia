%% function eroded_mask = erode3d(mask, noisemap)
%
% Description: Function to carry out 3D mask erosion based on a noise-based 
%              threshold.
% Input
% _____
%   mask              : logical n-d array masking the local field map
%   noisemap          : standard deviation of field map estimate (with the 
%                       same dimensions as the mask)
%
% Output
% ______
%   eroded_mask       : eroded mask (logical)
%
% Based on a design by Anita Karsa (UCL) Rewritten from an 
% implementation of Oliver Kiersnowski (UCL)
% Author: Patrick Fuchs @ UA
% patrick.fuchs@uantwerpen.be
% Date created: 23 July 2025
%
% see also SEPIAIOWRAPPER

function eroded_mask = erode3d(mask, noise_map)

arguments
    mask logical
    noise_map double
end

eroded_mask = zeros(size(mask));

noise_mask = 1./noise_map;
noise_mask(isnan(noise_mask) | isinf(noise_mask)) = 0;

thresh     = mean(noise_mask(noise_mask~=0)); 

eroded_mask(noise_mask >= thresh) = 1;
eroded_mask = eroded_mask && mask;

end
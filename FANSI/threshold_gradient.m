% Generate a binary image that acts as regularization weight, based on the gradient
% of the magnitude data.
%
% Based on the gradient_mask function in MEDI.
% Created by Ildar Khalidov in 2010
%
% Last modified by Carlos Milovic in 2017.03.30
%

function [ wG ] = threshold_gradient( gm, mask, noise, percentage )
%
% input:
% gm - gradient image (of the magnitude image, or T2* map, etc), as vector or magnitude.
% mask - binary 3D image that defines the ROI.
% noise - estimaded noise standard deviation in the complex signal
%
% output:
% wG - binary image. 0 denotes a relevant gradient, 1 elsewhere.
%

field_noise_level = noise;
denominator = sum(mask(:)==1);
wG = gm;
numerator = sum(wG(:)>field_noise_level);
if  (numerator/denominator)>percentage
    while (numerator/denominator)>percentage
        field_noise_level = field_noise_level*1.05;
        numerator = sum(wG(:)>field_noise_level);
    end
else
    while (numerator/denominator)<percentage
        field_noise_level = field_noise_level*.95;
        numerator = sum(wG(:)>field_noise_level);
    end
end

wG = (wG<=field_noise_level);



end


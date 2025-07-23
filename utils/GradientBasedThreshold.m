%% function MaskThresholded = GradientBasedThreshold(FieldMap, Mask, Lambda)
%
% Description: Compute a mask based on the (magnitude) of the gradient of 
%              the fieldmap.
% Input
% _____
%   FieldMap          : local field map in [ppm] (3D)
%   Mask              : input brain mask (logical, dimensions as FieldMap)
%   Lambda            : threshold at which to remove the voxels from the
%                       input mask
%
% Output
% ______
%   MaskThresholded   : eroded mask (logical)
%
% Created by: Oliver C. Kiersnowski @ UCL
% Date created: 1 February 2023
% Modified for SEPIA by: Patrick Fuchs @ UA
% Date modified: 23 July 2025
    

function MaskThresholded = GradientBasedThreshold(FieldMap, Mask, Lambda)
    
    [Gx,Gy,Gz] = imgradientxyz(FieldMap);
    Gmag = sqrt(Gx.^2 + Gy.^2 + Gz.^2);
    
    MeanGmag = mean(Gmag(:));
    StdGmag = std(Gmag(:));
    
    Thresh = MeanGmag + Lambda*StdGmag;
    
    MaskThresholded = Mask;
    MaskThresholded(Gmag > Thresh) = 0;

end
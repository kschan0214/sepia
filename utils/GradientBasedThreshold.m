%% function [MaskThresholded, Gmag, Stats] = GradientBasedThreshold(FieldMap, Mask, Lambda)
%
% Description: Compute a mask based on the (magnitude) of the gradient of
%              the fieldmap.
% Input
% _____
%   FieldMap          : local field map in [Hz] (3D)
%   Mask              : input brain mask (logical, dimensions as FieldMap)
%   Lambda            : threshold at which to remove the voxels from the
%                       input mask
%
% Output
% ______
%   MaskThresholded   : eroded mask (logical)
%   Gmag              : magnitude of the field map gradient (same size as
%                       FieldMap), returned so callers can inspect/save it
%                       alongside the mask
%   Stats             : struct with fields mean, std, lambda, threshold
%                       (all in Hz/voxel except lambda) - the values
%                       actually used to compute the threshold, so callers
%                       can record them for traceability (e.g. in the
%                       refined mask's JSON sidecar)
%
% Created by: Oliver C. Kiersnowski @ UCL
% Date created: 1 February 2023
% Modified for SEPIA by: Patrick Fuchs @ UA
% Date modified: 23 July 2025


function [MaskThresholded, Gmag, Stats] = GradientBasedThreshold(FieldMap, Mask, Lambda)

    [Gx,Gy,Gz] = imgradientxyz(FieldMap);
    Gmag = sqrt(Gx.^2 + Gy.^2 + Gz.^2);

    MeanGmag = mean(Gmag(:));
    StdGmag = std(Gmag(:));

    Thresh = MeanGmag + Lambda*StdGmag;

    fprintf('Gradient magnitude of the field map: mean = %g Hz/voxel, std = %g Hz/voxel\n', MeanGmag, StdGmag);
    fprintf('Lambda = %g, threshold = mean + lambda*std = %g Hz/voxel\n', Lambda, Thresh);

    Stats = struct('mean', MeanGmag, 'std', StdGmag, 'lambda', Lambda, 'threshold', Thresh);

    MaskThresholded = Mask;
    MaskThresholded(Gmag > Thresh) = 0;

end
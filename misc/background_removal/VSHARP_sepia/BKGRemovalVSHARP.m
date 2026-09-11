%% function [RDF,mask]=BKGRemovalVSHARP(totalField,mask,matrixSize,voxelSize,varargin)
%
% Input
% _____
%   totalField      : total field
%   mask            : ROI mask
%   matrixSize      : image matrix size
%   voxelSize       : voxel dimensions, in mm (1x3 vector); default = [1 1 1],
%                     i.e. 'radius' is treated as being in voxel unit
%   varargin        : flags with
%       'radius'    -   vector of radii being used, in mm
%
% Ouput
% _____
%   RDF             : local field
%
% Description: compute QSM based on iterative LSQR
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 19 July 2017
% Date modified: 1 April 2019
% Date modified: 12 April 2019
% Date modified: 23 August 2026 (bug fix: added missing k-space deconvolution step,
%                                which was causing incomplete background field removal
%                                radius is now specified in mm instead of voxel;
%                                added voxelSize input to support anisotropic voxels)
%
function [RDF,mask]=BKGRemovalVSHARP(totalField,mask,matrixSize,voxelSize,varargin)
if nargin < 4 || isempty(voxelSize)
    voxelSize = [1 1 1];
end
% parse argument input
[radius, threshold] = parse_varargin_VSHARP(varargin);

% zero padding
totalField  = padarray(totalField,[1 1 1],'both');
mask        = padarray(mask,[1 1 1],'both');
matrixSize  = size(totalField);

% total field in k-space
kTotalField = fftn(totalField);

DiffMask = zeros([matrixSize, length(radius)], 'like', totalField);
Mask_Sharp = zeros([matrixSize, length(radius)], 'like', totalField);
Del_Sharp = zeros([matrixSize, length(radius)], 'like', totalField);
RDF = 0;
% variable kernel size
for k = 1:length(radius)
    % get radius, in mm
    radiusCurrent = radius(k);

    % Sphere kernel in k-space (radius in mm, converted internally per axis using voxelSize)
    sphereKernel = SphereKernel(matrixSize,radiusCurrent,voxelSize);

    % erode mask to remove convolution artifacts
    % number of voxels along each axis needed to contain a sphere of
    % radiusCurrent mm, matching SphereKernel's kernel window
    erode_size = 2*ceil(radiusCurrent ./ voxelSize(:).') + 1;
    msk_sharp = imerode(mask, strel('line', erode_size(2), 0));
    msk_sharp = imerode(msk_sharp, strel('line', erode_size(1), 90));
    msk_sharp = permute(msk_sharp, [1,3,2]);
    msk_sharp = imerode(msk_sharp, strel('line', erode_size(3), 0));
    msk_sharp = permute(msk_sharp, [1,3,2]);

    Mask_Sharp(:,:,:,k) = msk_sharp; 
    Del_Sharp(:,:,:,k) = sphereKernel; 
    
    if k == 1
        DiffMask(:,:,:,1) = Mask_Sharp(:,:,:,1);
    else
        % boundary voxels between two spheres
        DiffMask(:,:,:,k) = Mask_Sharp(:,:,:,k) - Mask_Sharp(:,:,:,k-1);
    end
    % forward SHARP operator: subtract the spherical mean value
    fieldSharp = ifftn(Del_Sharp(:,:,:,k) .* kTotalField);

    % mask out unreliable boundary data before deconvolution (as in SHARP.m)
    fieldSharp = fieldSharp .* msk_sharp;

    % deconvolution in k-space to recover the true local field for this kernel size
    kernelCurrent   = Del_Sharp(:,:,:,k);
    kRDF            = fftn(fieldSharp) ./ kernelCurrent;
    kRDF(abs(kernelCurrent) <= threshold) = 0;
    kRDF(isnan(kRDF)) = 0;
    kRDF(isinf(kRDF)) = 0;

    % if k~=1, put back the phase in the boundary between current kernel
    % and previous kernel
    RDF = RDF + DiffMask(:,:,:,k) .* real(ifftn(kRDF));
end
%  largest mask
mask = Mask_Sharp(:,:,:,end);     

RDF = RDF .* mask ;

RDF = RDF(2:end-1,2:end-1,2:end-1);
mask = mask(2:end-1,2:end-1,2:end-1);

end

%% parser
function [radius, threshold] = parse_varargin_VSHARP(arg)
radius    = 5:-1:1;
threshold = 0.05;
if ~isempty(arg)
    for kvar = 1:length(arg)
        if strcmpi(arg{kvar},'radius')
            tmp = arg{kvar+1};
            radius = sort(tmp,'descend');
        end
        if strcmpi(arg{kvar},'threshold')
            threshold = arg{kvar+1};
        end
    end
end
end
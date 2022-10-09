% Sophisticated Harmonic Artifact Reduction for Phase data (SHARP)
%   RDF = SHARP(iFreq, Mask,matrix_size,voxel_size, radius,threshold)
% 
%   output
%   RDF - the relative difference field, or local field
% 
%   input
%   iFreq - the unwrapped field map
%   Mask - a binary 3D matrix denoting the Region Of Interest
%   matrix_size - the size of the 3D matrix
%   voxel_size - the size of the voxel in mm
%   radius (optional) - the radius of the spherical mean value operation
%   treshodl (optional) - the threshold used in Truncated SVD
%
%   When using the code, please cite 
%   F. Schweser et al. NeuroImage 54 (2011)2789?2807
%
%   Created by Tian Liu in 2010
%   Last modified by Tian Liu on 2011.02.01

function RDF=SHARP(iFreq,Mask, matrix_size,voxel_size,radius,threshold)

if (nargin<6)
    threshold=0.00;
end
if (nargin<5)
    radius=round(6/max(voxel_size)) * max(voxel_size);
end


%erode the Mask to exclude the boundary
Mask = SMV(Mask, matrix_size,voxel_size,radius)>0.999; 

% generate the convolution/deconvolution kernel
kernel = SMV_kernel(matrix_size, voxel_size, radius);

% subtraction of the spherical mean value
iFreq = ifftn(fftn(iFreq).*kernel);

% mask out the unreliable phase data
iFreq_mask = iFreq.*Mask;


% deconvolution in k-space
kRDF = fftn(iFreq_mask)./kernel;
% truncation in k-space
kRDF(abs(kernel)<=threshold) = 0;
kRDF(isnan(kRDF)) = 0;
kRDF(isinf(kRDF)) = 0;
% results in image space
RDF = ifftn(kRDF).*Mask;


% Regularization Enabled Sophisticated Harmonic Artifact Reduction for Phase data (RESHARP)
%   RDF = RESHARP(iFreq, Mask,matrix_size,voxel_size, radius,alpha)
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
%   alpha (optional) - the regularizaiton parameter used in Tikhonov
%
%   When using the code, please cite 
%   Sun et al. MRM 2014;71(3):1151-7
%   Created by Tian Liu on 2014.02.01


function RDF=RESHARP(iFreq,Mask, matrix_size,voxel_size,radius,alpha)

if (nargin<6)
    alpha=0.01;
end
if (nargin<5)
    radius=round(6/max(voxel_size)) * max(voxel_size);
end

% generate the convolution/deconvolution kernel
S = SMV_kernel(matrix_size, voxel_size, radius);
M1 = SMV(Mask, matrix_size, voxel_size, radius)>0.999;

A = @(x) M1.*ifftn(S.*fftn(x));
Ah = @(x) ifftn(S.*fftn(M1.*x));
b = Ah(A(iFreq));
FW = @(x) (A(Ah(x))+alpha*(x));

RDF = cgsolve(FW, b, 1e-2,40,0);


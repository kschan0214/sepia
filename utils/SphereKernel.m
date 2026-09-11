%% function sphereKernel = SphereKernel(matrixSize,radius,voxelSize)
%
% Input
% --------------
%   matrixSize      : chi matrix size
%   radius          : radius of sphere, in mm
%   voxelSize       : voxel dimensions, in mm (1x3 vector); default = [1 1 1],
%                     i.e. radius is treated as being in voxel unit
%
% Output
% --------------
%   sphereKernel    : sphere in k-space with matrix size matched with chi
%
% Description: Sphere kernel in k-space, modified based on Bilgic's code.
%              The sphere is defined in physical (mm) space so that, with
%              anisotropic voxel spacing, the kernel is an ellipsoid in
%              voxel-index space that corresponds to a true sphere of the
%              specified radius in physical units.
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 19 July 2017
% Date modified: 23 August 2026 (radius now in mm; supports anisotropic voxelSize)
%
%
function sphereKernel = SphereKernel(matrixSize,radius,voxelSize)

if nargin < 3
    voxelSize = [1 1 1];
end

% number of voxels needed along each axis to fully contain the sphere
halfwidth = ceil(radius ./ voxelSize(:).');

% create a 3D grid in physical (mm) space, centred at the origin
% (dimension order follows matrixSize: dim1-dim2-dim3)
[Y,X,Z] = meshgrid(-halfwidth(2):halfwidth(2), -halfwidth(1):halfwidth(1), -halfwidth(3):halfwidth(3));
X = X * voxelSize(1);
Y = Y * voxelSize(2);
Z = Z * voxelSize(3);

% modify 3D grid to sphere (in physical space)
sphere = (X.^2 + Y.^2 + Z.^2) <= radius^2;

% normalised the kernel values such that the sum of the kernel = 1
sphere = -sphere / sum(sphere(:));
centreIdx = halfwidth + 1;
sphere(centreIdx(1),centreIdx(2),centreIdx(3)) = 1 + sphere(centreIdx(1),centreIdx(2),centreIdx(3));

% match the matrix size of the spherical kernel to the image size
Kernel = zeros(matrixSize);
Kernel( 1+matrixSize(1)/2 - halfwidth(1) : 1+matrixSize(1)/2 + halfwidth(1), ...
        1+matrixSize(2)/2 - halfwidth(2) : 1+matrixSize(2)/2 + halfwidth(2), ...
        1+matrixSize(3)/2 - halfwidth(3) : 1+matrixSize(3)/2 + halfwidth(3) ) = sphere;

% Spherical kernel in k-space
sphereKernel = fftn(fftshift(Kernel));

end

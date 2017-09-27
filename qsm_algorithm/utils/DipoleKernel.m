%% function [dipoleKernel,dKComponents] = DipoleKernel(matrixSize,voxelSize,b0dir)
%
% Description: Create dipole kernel in k-space with input matrix dimensions
%              and spatial resolution
% Input
% _____
%   matrixSize        : image matrix size
%   voxelSize         : spatial resolution of image 
%   b0dir             : static magnetic field direction (optional)
%
% Output
% ______
%   dipoleKernal      : Dipole kernel (in k-space)
%   dKComponents      : dipole kernel components
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 24 March 2017
% Date last modified: 27 September 2017
%
function [dipoleKernel,dKComponents] = DipoleKernel(matrixSize,voxelSize,b0dir)

if isempty(b0dir)
    b0dir = [0 0 1];
end

% KC: create 3D matrix in k-space
% [ky,kx,kz] = ndgrid(-matrixSize(1)/2:matrixSize(1)/2-1, ...
%                     -matrixSize(2)/2:matrixSize(2)/2-1, ...
%                     -matrixSize(3)/2:matrixSize(3)/2-1);
[ky,kx,kz] = meshgrid(-matrixSize(2)/2:matrixSize(2)/2-1, ...
                      -matrixSize(1)/2:matrixSize(1)/2-1, ...
                      -matrixSize(3)/2:matrixSize(3)/2-1);

% KC: assign k-vectors
kx = (kx / max(abs(kx(:)))) / voxelSize(1);
ky = (ky / max(abs(ky(:)))) / voxelSize(2);
kz = (kz / max(abs(kz(:)))) / voxelSize(3);

k2 = kx.^2 + ky.^2 + kz.^2;

% KC: shift the centre of k-space to matrix corners
% KC: second term represents (cos beta).^2 where beta is the angle between
%     k-vector and the static B field
% dipoleKernel = fftshift( 1/3 - (kz ).^2 ./ (k2 + eps) );
% 260917: correct for B0 direction, b0dir = [x,y,z]; 
% dipoleKernel = fftshift( 1/3 - (kx*b0dir(2) + ky*b0dir(1) + kz*b0dir(3)).^2 ./ (k2 + eps) );
dipoleKernel = fftshift( 1/3 - (kx*b0dir(1) + ky*b0dir(2) + kz*b0dir(3)).^2 ./ (k2 + eps) );

dKComponents.kx = kx;
dKComponents.ky = ky;
dKComponents.kz = kz;
dKComponents.k2 = k2;

end
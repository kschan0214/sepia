%% function dipoleKernal = DipoleKernal(matrixSize,voxelSize)
%
% Description: Create dipole kernel in k-space with input matrix dimensions
%              and spatial resolution
% Input
% _____
%   matrixSize        : image matrix size
%   voxelSize         : spatial resolution of image 
%
% Output
% ______
%   dipoleKernal                : Dipole kernel (in k-space)
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 24 March 2017
% Date last modified: 28 June 2017
%
function dipoleKernal = DipoleKernal(matrixSize,voxelSize)

% KC: create 3D matrix in k-space
[kx,ky,kz] = ndgrid(-matrixSize(1)/2:matrixSize(1)/2-1, ...
                    -matrixSize(2)/2:matrixSize(2)/2-1, ...
                    -matrixSize(3)/2:matrixSize(3)/2-1);

% KC: assign k-vectors
kx = (kx / max(abs(kx(:)))) / voxelSize(1);
ky = (ky / max(abs(ky(:)))) / voxelSize(2);
kz = (kz / max(abs(kz(:)))) / voxelSize(3);

k2 = kx.^2 + ky.^2 + kz.^2;

% KC: shift the centre of k-space to matrix corners
% KC: second term represents (cos beta).^2 where beta is the angle between
%     k-vector and the static B field
dipoleKernal = fftshift( 1/3 - (kz ).^2 ./ (k2 + eps) );

end
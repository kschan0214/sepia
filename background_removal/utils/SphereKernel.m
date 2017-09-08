%% function sphereKernel = SphereKernel(matrixSize,radius)
%
% Input
% --------------
%   matrixSize      : chi matrix size
%   radius          : radisu of sphere (in order of pixel)
%
% Output
% --------------
%   sphereKernel    : sphere in k-space with matrix size matched with chi 
%
% Description: Sphere kernel in k-space, modified based on Bilgic's code
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 19 July 2017
% Date last modified:
%
%
function sphereKernel = SphereKernel(matrixSize,radius)

% create 3D square    
[a,b,c] = meshgrid(-radius:radius, -radius:radius, -radius:radius);

% modify 3D square to sphere
sphere = (a.^2 / radius^2 + b.^2 / radius^2 + c.^2 / radius^2 ) <= 1;

% normalised the kernel values such that the sum of the kernel = 1
sphere = -sphere / sum(sphere(:));
sphere(radius+1,radius+1,radius+1) = 1 + sphere(radius+1,radius+1,radius+1);

% match the matrix size of the shperical kernel to the image size
Kernel = zeros(matrixSize);
Kernel( 1+matrixSize(1)/2 - radius : 1+matrixSize(1)/2 + radius, ...
        1+matrixSize(2)/2 - radius : 1+matrixSize(2)/2 + radius, ...
        1+matrixSize(3)/2 - radius : 1+matrixSize(3)/2 + radius ) = sphere;

% Spherical kernel in k-space
sphereKernel = fftn(fftshift(Kernel));    

end
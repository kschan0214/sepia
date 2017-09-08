% Generate a Spherical kernel with the sum normalized to one
%   y = SMV_kernel(matrix_size,voxel_size, radius)
%   
%   output
%   y - kernel
% 
%   input
%   matrix_size - the dimension of the field of view
%   voxel_size - the size of the voxel in mm
%   radius - the raidus of the sphere in mm
%
%   Created by Tian Liu in 2010
%   Modified by Tian on 2011.02.01
%   Modified by Tian on 2011.03.14 The sphere is now rendered.
%   Last modified by Tian Liu on 2013.07.23

function y = sphere_kernel(matrix_size,voxel_size, radius)

[Y,X,Z]=meshgrid(-matrix_size(2)/2:matrix_size(2)/2-1,...
                 -matrix_size(1)/2:matrix_size(1)/2-1,...
                 -matrix_size(3)/2:matrix_size(3)/2-1);

X = X*voxel_size(1);
Y = Y*voxel_size(2);
Z = Z*voxel_size(3);
Sphere_out = (   max(abs(X)-0.5*voxel_size(1),0).^2 ... 
               +max(abs(Y)-0.5*voxel_size(2),0).^2 ...
               +max(abs(Z)-0.5*voxel_size(3),0).^2 )>radius^2;

Sphere_in = (  (abs(X)+0.5*voxel_size(1)).^2 ... 
                +(abs(Y)+0.5*voxel_size(2)).^2 ...
                +(abs(Z)+0.5*voxel_size(3)).^2 )<=radius^2; 

            
Sphere_mid = zeros(matrix_size);

split = 10; %such that error is controlled at <1/(2*10)
[X_v Y_v Z_v] = meshgrid(-split+0.5:split-0.5, -split+0.5:split-0.5, -split+0.5:split-0.5);
X_v = X_v/(2*split);
Y_v = Y_v/(2*split);
Z_v = Z_v/(2*split);

shell = 1-Sphere_in-Sphere_out;
X = X(shell==1);
Y = Y(shell==1);
Z = Z(shell==1);
shell_val = zeros(size(X));

for i = 1:length(X)
    xx = X(i);
    yy = Y(i);
    zz = Z(i);
   
    occupied = ( (xx+X_v*voxel_size(1)).^2+...
        (yy+Y_v*voxel_size(2)).^2+...
        (zz+Z_v*voxel_size(3)).^2)<=radius^2;
    shell_val(i) = sum(occupied(:))/numel(X_v);
end

Sphere_mid(shell==1) = shell_val;

Sphere = Sphere_in+Sphere_mid;    
Sphere = Sphere/sum(Sphere(:));
y = fftn(fftshift(Sphere));
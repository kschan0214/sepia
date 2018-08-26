% Generate a kernel that performs the removal of the spherical mean value
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


function y = SMV_kernel(matrix_size,voxel_size, radius)

y = 1-sphere_kernel(matrix_size, voxel_size,radius);
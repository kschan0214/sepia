% Unwrap phase using Laplacian operation
% Schofield and Zhu, 2003 M.A. Schofield and Y. Zhu, Fast phase unwrapping
% algorithm for interferometric applications, Opt. Lett. 28 (2003), pp. 1194?196
%   [iFreq ] = unwrapLaplacian(iFreq_raw, matrix_size, voxel_size)
% 
%   output
%   iFreq - unwrapped phase
%
%   input
%   iFreq_raw - A wrapped field map
%   matrix_size - dimension of the field of view
%   voxel_size - size of voxel in mm
%   
%   Created by Tian Liu on 2008.8.10
%   Modified by Tian Liu on 2011.01.26
%   Last modified by Tian Liu on 2013.07.24
%   Last modified by Tian Liu on 2014.02.10


function [iFreq ] = unwrapLaplacian(iFreq_raw, matrix_size, voxel_size)


if (nargin<3)
    voxel_size = [1 1 1];
end

if length(matrix_size)==2
    matrix_size(3) = 1;
end

if length(voxel_size)==2
    voxel_size(3) = 1;
end


[Y,X,Z]=meshgrid(-matrix_size(2)/2:matrix_size(2)/2-1,...
                 -matrix_size(1)/2:matrix_size(1)/2-1,...
                 -matrix_size(3)/2:matrix_size(3)/2-1);

X = X*voxel_size(1);
Y = Y*voxel_size(2);
Z = Z*voxel_size(3);


if matrix_size(3)>1
    h=( (X==0) & (Y==0)&(Z==0) ) .*1;
    k=6 * del2(h,voxel_size(1), voxel_size(2), voxel_size(3));
    kernel=fftn(fftshift(k));
else 
    h=( (X==0) & (Y==0) ) .*1;
    k=4* del2(h,voxel_size(1), voxel_size(2));
    kernel=fftn(fftshift(k));
end

inv_kernel = 1./kernel;
inv_kernel(isinf(inv_kernel))=0;
inv_kernel(abs(kernel)<1e-10)=0;



first_term  = cos(iFreq_raw) .*  ifftn(  kernel.* fftn(sin(iFreq_raw)));  %  
second_term = sin(iFreq_raw) .*  ifftn( kernel.* fftn(cos(iFreq_raw)));  %  
phi_est   =  ifftn( inv_kernel.* fftn(  (first_term - second_term) ));   % 
iFreq = phi_est;




end




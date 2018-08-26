% discrete gradient using central difference
%   Gx = cdiv(x, voxel_size)
%   
%   output
%   Gx - the gradient
% 
%   input
%   x - the scalar field
%   voxel_size - the size of the voxel
%
%   Created by Ludovic de Rochefort in 2009
%   Changed by Shuai Wang and Tian Liu on 2011.03.28
%   Last modified by Tian Liu on 2013.07.24

function Gx = cgrad(x, voxel_size )
if(nargin <2)
    voxel_size = [1 1 1];
end

Dx = -0.5*circshift(x,[1 0 0]) + 0.5*circshift(x,[-1 0 0]);
Dx(1,:,:) = -x(1,:,:)+x(2,:,:);
Dx(end,:,:) = x(end,:,:) - x(end-1,:,:);
Dx = Dx/voxel_size(1);

Dy = -0.5*circshift(x,[0 1 0]) + 0.5*circshift(x,[0 -1 0]);
Dy(:,1,:) = -x(:,1,:)+x(:,2,:);
Dy(:,end,:) = x(:,end,:) - x(:,end-1,:);
Dy = Dy/voxel_size(2);

Dz = -0.5*circshift(x,[0 0 1]) + 0.5*circshift(x,[0 0 -1]);
Dz(:,:,1) = -x(:,:,1)+x(:,:,2);
Dz(:,:,end) = x(:,:,end) - x(:,:,end-1);
Dz = Dz/voxel_size(3);

Gx=cat(4,Dx,Dy,Dz);
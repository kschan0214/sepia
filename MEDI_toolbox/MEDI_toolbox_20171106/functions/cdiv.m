% discrete divergence using central difference
%   c = cdiv(Gx, voxel_size)
%   
%   output
%   c is the divergence
% 
%   input
%   Gx is a vector field
%   voxel_size is the size of the voxel
%
%   Created by Ludovic de Rochefort in 2009
%   Changed by Shuai Wang and Tian Liu on 2011.03.28
%   Last modified by Tian Liu on 2013.07.24

function c = cdiv(Gx, voxel_size)

cx = 0.5*circshift(Gx(:,:,:,1),[1 0 0]) - 0.5*circshift(Gx(:,:,:,1),[-1 0 0]);
cx(1,:,:) = -Gx(1,:,:,1)-0.5*Gx(2,:,:,1);
cx(2,:,:) = Gx(1,:,:,1)-0.5*Gx(3,:,:,1);
cx(end-1,:,:) = 0.5*Gx(end-2,:,:,1) - Gx(end,:,:,1);
cx(end,:,:) = 0.5*Gx(end-1,:,:,1) + Gx(end,:,:,1);
cx = cx/voxel_size(1);

cy = 0.5*circshift(Gx(:,:,:,2),[0 1 0]) - 0.5*circshift(Gx(:,:,:,2),[0 -1 0]);
cy(:,1,:) = -Gx(:,1,:,2)-0.5*Gx(:,2,:,2);
cy(:,2,:) = Gx(:,1,:,2)-0.5*Gx(:,3,:,2);
cy(:,end-1,:) = 0.5*Gx(:,end-2,:,2) - Gx(:,end,:,2);
cy(:,end,:) = 0.5*Gx(:,end-1,:,2) + Gx(:,end,:,2);
cy = cy/voxel_size(2);

cz = 0.5*circshift(Gx(:,:,:,3),[0 0 1]) - 0.5*circshift(Gx(:,:,:,3),[0 0 -1]);
cz(:,:,1) = -Gx(:,:,1,3)-0.5*Gx(:,:,2,3);
cz(:,:,2) = Gx(:,:,1,3)-0.5*Gx(:,:,3,3);
cz(:,:,end-1) = 0.5*Gx(:,:,end-2,3) - Gx(:,:,end,3);
cz(:,:,end) = 0.5*Gx(:,:,end-1,3) + Gx(:,:,end,3);
cz = cz/voxel_size(3);

c=cx+cy+cz;
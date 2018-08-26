function res = adjD(voxel_size,y)

res = (y([1,1:end-1],:,:,1) - y(:,:,:,1))/voxel_size(1) + ...
    (y(:,[1,1:end-1,],:,2) - y(:,:,:,2))/voxel_size(2) + ...
    (y(:,:,[1,1:end-1],3) - y(:,:,:,3))/voxel_size(3);

return;


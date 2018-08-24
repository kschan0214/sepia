function res = D(voxel_size,image)


res=zeros(size(image,1),size(image,2),size(image,3),3);
res(:,:,:,1) = (image([2:end,end],:,:) - image)/voxel_size(1);
res(:,:,:,2) = (image(:,[2:end,end],:) - image)/voxel_size(2);
res(:,:,:,3) = (image(:,:,[2:end,end]) - image)/voxel_size(3);


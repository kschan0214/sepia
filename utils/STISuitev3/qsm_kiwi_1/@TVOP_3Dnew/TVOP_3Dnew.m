function  res = TVOP_3Dnew(voxel_size,mask)

res.voxel_size = voxel_size;
res.mask = mask;
res.adjoint = 0;
res = class(res,'TVOP_3Dnew');


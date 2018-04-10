function res = mtimes(a,b)

if a.adjoint
	res = adjD(a.voxel_size,b.*a.mask);
else
	res = D(a.voxel_size,b).*a.mask;
end




    

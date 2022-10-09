% handy function to save variable in nifti format
function save_nii_quick(templateHeader,img,filename,datatype)

if nargin < 4
    datatype = 16;
end

nii = templateHeader;
nii.img = single(img);
nii.hdr.dime.datatype = datatype;
nii.hdr.dime.dim(5) = size(img,4);
nii.hdr.dime.dim(1) = ndims(img);

% assume the input image contains the true values
nii.hdr.dime.scl_inter = 0;
nii.hdr.dime.scl_slope = 1;

save_untouch_nii(nii,filename);

end
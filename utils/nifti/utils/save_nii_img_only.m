%% function save_nii_img_only(headerfilename,savefilename,images)
%
% Description: Wrapper script for saving with target image header
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 24 January 2017
% Date last modified:
%
function save_nii_img_only(headerfilename,savefilename,img,datatype)

if nargin < 4
    datatype = 16;
end

nii = load_untouch_nii(headerfilename);
nii.img = single(img);
nii.hdr.dime.datatype = datatype;
nii.hdr.dime.dim(5) = size(img,4);
nii.hdr.dime.pixdim(isnan(nii.hdr.dime.pixdim)) = 1;
if size(img,4) > 1
    nii.hdr.dime.dim(1) = 4;
end

% assume the input image contains the true values
nii.hdr.dime.scl_inter = 0;
nii.hdr.dime.scl_slope = 1;

% if nii.hdr.dime.datatype ==4
%     nii.hdr.dime.datatype = 16;
% end
% nii.img = images;

save_untouch_nii(nii,savefilename);
end
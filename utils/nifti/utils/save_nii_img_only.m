%% function save_nii_img_only(headerfilename,savefilename,images)
%
% Description: Wrapper script for saving with target image header
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 24 January 2017
% Date last modified:
%
function save_nii_img_only(headerfilename,savefilename,images)
load_module_NIfTI;
nii = load_untouch_nii(headerfilename);
nii.img = images;
save_untouch_nii(nii,savefilename);
end
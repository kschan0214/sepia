%% function img = load_nii_img_only(filename)
%
% Description: Load images data only
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 9 October 2016
% Date last modified:
%
function img = load_nii_img_only(filename)
load_module_NIfTI;
a = load_untouch_nii(filename);
img = a.img;
end
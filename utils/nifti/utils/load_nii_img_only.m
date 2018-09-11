%% function img = load_nii_img_only(filename)
%
% Description: Load images data only
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 9 October 2016
% Date last modified: 11 September 2018
%
function img = load_nii_img_only(filename)
a = load_untouch_nii(filename);
img = a.img;
end
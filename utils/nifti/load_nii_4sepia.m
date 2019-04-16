%% function output = function_name(input)
%
% Usage:
%
% Input
% --------------
%
% Output
% --------------
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 16 April 2019
% Date last modified:
%
%
function nii = load_nii_4sepia(filename, img_idx, dim5_idx, dim6_idx, dim7_idx, ...
			old_RGB, slice_idx)

nii = load_untouch_nii(filename, img_idx, dim5_idx, dim6_idx, dim7_idx, ...
			old_RGB, slice_idx);

end
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
try 
    info = niftiinfo(filename);
    img = double(niftiread(info))*info.MultiplicativeScaling + info.AdditiveOffset;
catch
    a = load_untouch_nii(filename);
    img = double(a.img)*a.hdr.dime.scl_slope + a.hdr.dime.scl_inter;
end
end
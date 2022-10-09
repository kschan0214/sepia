%% function img = load_nii_img_only(filename)
%
% Description: Load images data only
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 9 October 2016
% Date modified: 11 September 2018
% Date modified: 26 June 2020 (v0.8.0)
%
function img = load_nii_img_only(filename)
% try 
%     info = niftiinfo(filename);
%     if info.MultiplicativeScaling == 0
%         img = double(niftiread(info));
%     else
%         img = double(niftiread(info))*info.MultiplicativeScaling + info.AdditiveOffset;
%     end
% catch
    a = load_untouch_nii(filename);
    if a.hdr.dime.scl_slope == 0
        img = double(a.img);
    else
        img = double(a.img)*a.hdr.dime.scl_slope + a.hdr.dime.scl_inter;
    end
% end
end
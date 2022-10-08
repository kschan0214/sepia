%% img_hybrid = compute_hybrid_t1w_chi(t1w, chi, mask, mu, thres )
%
% Input
% --------------
% t1w           : T1w image (or R1)
% chi           : Chimap in T1w space, in ppm
% mask          : brain mask;
%
% Output
% --------------
% img_hybrid    : hybrid contrast 
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 7 October 2022
% Date last modified:
%
%
function img_hybrid = compute_hybrid_t1w_chi(t1w, chi, mask, mu, thres )

if nargin < 5
    thres = 0.5; % percent
end
if nargin < 4
    mu = 400; % a.u.
end

t1w_norm = (t1w - prctile(t1w(mask>0),thres)) / (prctile(t1w(mask>0),100 - thres) - prctile(t1w(mask>0),thres));
t1w_norm(t1w_norm>1) = 1; t1w_norm(t1w_norm<0) = 0; t1w_norm = t1w_norm*255;

img_hybrid = t1w_norm - mu*chi;
img_hybrid(img_hybrid<0) = 0; img_hybrid(img_hybrid>255) = 255;

end
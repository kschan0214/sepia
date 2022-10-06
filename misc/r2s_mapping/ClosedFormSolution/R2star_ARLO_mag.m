%% [r2s, t2s, m0] = R2star_ARLO_mag(img,te,m0mode)
%
% Input
% --------------
% img           : multiecho images, time in last dimension
% te            : echo times
% m0mode        : method to compute m0, 'default','weighted','average'
%
% Output
% --------------
%   r2s         : R2* map
%   t2s         : T2* map
%   m0          : M0, T1-weighted image
%
% Description: Fast monoexponential fitting of R2* by auto-regression
% 
% ref: Pei et al. MRM 73:843-850(2015)
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 7 October, 2016
% Date last modified: 29 June 2017
%
function [r2s, t2s, s0] = R2star_ARLO_mag(img,te,s0mode)

% set range of R2* and T2*
minT2s = min(te)/20;
maxT2s = max(te)*20;
ranget2s = [minT2s, maxT2s];
ranger2s = [1/maxT2s, 1/minT2s];

% set m0 extrapolation method
if nargin < 3
    s0mode = '1st echo';
end

% disgard phase information
img = double(abs(img));
te = double(te);


%% main
% compute T2*  by ARLO
t2s = ARLO(img,te);
% convert T2* to R2*
r2s = 1./t2s;

% set range
t2s = SetImgRange(t2s,ranget2s);
r2s = SetImgRange(r2s,ranger2s);

% calculate s0
s0 = ComputeM0GivenR2star(r2s,te,img,s0mode);

end

%% uility function
function res = SetImgRange(img,range)
    imgMax = range(2);
    imgMin = range(1);
    img(img<imgMin) = imgMin;
    img(img>imgMax) = imgMax;
    img(isinf(img)) = imgMin;
    img(isnan(img)) = imgMin;
    res = img;
end
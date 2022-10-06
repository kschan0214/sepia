%% [r2s,t2s,m0] = R2star_GS(img,te,m0mode)
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
% Description: Monoexponential R2* computation by geometric sum
%              Only works for even echo-spacing data
% Details check notebook
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 5 February 2017
% Date last modified: 30 June 2017
%
function [r2s,t2s,s0] = R2star_GS(img,te,s0mode)

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

% get dimension of data
ndim = ndims(img);

% Get echo spacing
dt = te(2)-te(1);

%% Main
% Q = sum of all echoes, assuming time in 4th dimension
Q = sum(img,ndim);

% Get 1st echo image and last echo image
I1 = img(:,:,:,1);
Iend = img(:,:,:,end);

% Results of geometric product
r2s = real(log((Q-Iend)./(Q-I1))./ dt);

% convert to T2*
t2s = 1./r2s;

%% Result preparation
% set range
r2s = SetImgRange(r2s,ranger2s);
t2s = SetImgRange(t2s,ranget2s);

% calculate m0
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

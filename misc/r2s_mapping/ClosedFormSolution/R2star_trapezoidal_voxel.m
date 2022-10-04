%% [r2s,t2s,m0] = R2star_trapezoidal_voxel(s,te,m0mode)
%
% Input
% --------------
% img           : multiecho images, time in last dimension
% te            : echo times
% m0mode        : method to compute m0, 'default','weighted','average'
%
% Output
% --------------
%   r2s         : R2*
%   t2s         : T2*
%   m0          : M0, T1-weighted images
%
% Description: Trapezoidal approximation of integration
% Code From Jose Marques
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 13 October 2016
% Date last modified: 24 September 2017
%
function [r2s,t2s,s0] = R2star_trapezoidal_voxel(s,te,s0mode)

% set range of R2* and T2*
minT2s = min(te)/20;
maxT2s = max(te)*20;
ranget2s = [minT2s, maxT2s];
ranger2s = [1/maxT2s, 1/minT2s];

% set m0 extrapolation method
if nargin < 3
    s0mode = '1stecho';
end

% disgard phase information
s = abs(s);

%% main
temp=0;
for k=1:length(s)-1
    temp=temp+0.5*(s(k)+s(k+1))*(te(k+1)-te(k));
end
    % very fast estimation
t2s=temp./(s(1)-s(end));
    
r2s = 1./t2s;

%% Result preparation
% set range
r2s = SetImgRange(r2s,ranger2s);
t2s = SetImgRange(t2s,ranget2s);
    
% calculate m0
s0 = ComputeM0GivenR2star(r2s,te,s,s0mode);
    
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
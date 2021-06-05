%% [r2s,t2s,m0] = R2star_trapezoidal(img,te,m0mode)
%
% Input
% --------------
% img           : multiecho images, time in last dimension
% te            : echo times (in s)
% m0mode        : method to compute m0, 'default','weighted','average'
%
% Output
% --------------
%   r2s         : R2* map
%   t2s         : T2* map
%   m0          : M0, T1-weighted images
%
% Description: Trapezoidal approximation of integration
% Code From Jose P. Marques
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 13 October 2016
% Date modified: 24 September 2017
% Date last modified: 21 April 2018
%
function [r2s,t2s,s0] = R2star_trapezoidal(img,te,s0mode)

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
img = double(abs(img));
te = double(te);

dims = size(img);

%% main
% Trapezoidal approximation of integration
temp=0;
for k=1:dims(4)-1
    temp=temp+0.5*(img(:,:,:,k)+img(:,:,:,k+1))*(te(k+1)-te(k));
end

% very fast estimation
t2s=temp./(img(:,:,:,1)-img(:,:,:,end));
    
r2s = 1./t2s;

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
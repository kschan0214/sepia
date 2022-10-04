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
% Date created: 
% Date last modified:
%
%
function [r2s,t2s,m0] = R2starmapping_CRSI(img,te,m0mode,sigma,m)

% set range of R2* and T2*
minT2s = min(te)/20;
maxT2s = max(te)*20;
ranget2s = [minT2s, maxT2s];
ranger2s = [1/maxT2s, 1/minT2s];

% set m0 extrapolation method
if nargin < 3
    m0mode = 'default';
end

% p = abs(img).^2;
p = real(img).^2 + imag(img).^2;
nTE = length(te);

A = 0;
for n=1:nTE-1
    tmp=0;
    for km=0:m
        tmp = tmp + p(:,:,:,n).^((m-km+1)/(m+1)) .* p(:,:,:,n+1).^(km/(m+1));
    end
    A = A+((te(n+1)-te(n))/(m-1) * tmp);
end

r2s = (p(:,:,:,1) - p(:,:,:,end))./ (2*(A-2*sigma.^2*(te(end)-te(1))));
t2s = 1./r2s;
%% Result preparation
% set range
r2s = SetImgRange(r2s,ranger2s);
t2s = SetImgRange(t2s,ranget2s);
    
% calculate m0
m0 = ComputeM0GivenR2star(r2s,te,img,m0mode);

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

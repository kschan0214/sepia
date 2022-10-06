%% [r2s, t2s, m0] = R2star_pi(img,te,m0mode,method)
%
% Input
% --------------
% img           : multiecho images, time in last dimension
% te            : echo times
% m0mode        : method to compute m0, 'default','weighted','average'
% method        : method to compute R2* ('interleaved; or '1stecho')
%
% Output
% --------------
%   r2s         : R2*
%   t2s         : T2*
%   m0          : M0, T1-weighted images
%
% Description: R2* computation by using sequence of product
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 03 February 2017
% Date last modified: 30 June 2017
%
function [r2s, t2s, s0] = R2star_pi(img,te,s0mode,method)

% set range of R2* and T2*
minT2s = min(te)/20;
maxT2s = max(te)*20;
ranget2s = [minT2s, maxT2s];
ranger2s = [1/maxT2s, 1/minT2s];

% set m0 extrapolation method
if nargin < 3
    s0mode = '1st echo';
end

if nargin < 4
    method = 'interleaved';
end

% disgard phase information
img = double(abs(img));
te = double(te);

[nx, ny, nz, nt] = size(img);
%% Main
if strcmpi(method,'1st echo')
    tmp = img(:,:,:,2:end)./img(:,:,:,1);
	dt_sum = sum(te(2:end)-te(1));
else
    tmp = img(:,:,:,2:end)./img(:,:,:,1:end-1); 
	dt_sum = sum(te(2:end)-te(1:end-1));
end

Q = ones(nx,ny,nz);
for kt=1:nt-1
    Q = Q .* tmp(:,:,:,kt);
end 

r2s = -log(Q)/dt_sum;

%% Result preparation
% convert to T2*
t2s = 1./r2s;

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
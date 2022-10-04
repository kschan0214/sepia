%% [r2s,t2s,m0] = R2star_regression(img,te,mask,isParallel)
%
% Input
% --------------
% img           : multiecho images, time in last dimension
% te            : echo times
% mask          : signal mask
%
% Output
% --------------
% r2s           : R2* map
% t2s           : T2* map
% m0            : M0, T1-weighted image
%
% Description: R2* fitting by simple regression model
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 13 October 2017
% Date last modified:
%
%
function [r2s,t2s,m0] = R2star_regression(img,te,mask,isParallel)

% set range of R2* and T2*
minT2s = min(te)/20;
maxT2s = max(te)*20;
ranget2s = [minT2s, maxT2s];
ranger2s = [1/maxT2s, 1/minT2s];

if nargin < 4
    isParallel = false;
end

% if no mask input use default mask to speed up regression
if isempty(mask)
    mask = max(img,[],4)/(max(img(:))) >0.02;
end

img = double(img);
te = double(te);

[nx,ny,nz,nt] = size(img);

x = ones(nt,2);
x(:,2) = -te(:);
y = permute(log(abs(img)),[4 1 2 3]);

y = reshape(y,[size(y,1) numel(y)/size(y,1)]);
b = x\y;
r2s = reshape( b(2,:),nx,ny,nz) .* mask;
m0 = exp(reshape( b(1,:),nx,ny,nz)) .* mask;

% r2s = zeros(nx,ny,nz);
% m0 = zeros(nx,ny,nz);
% if isParallel
%     for kz=1:nz
%         for ky=1:ny
%             parfor kx=1:nx
%                 if mask(kx,ky,kz)==1
%                     b = x\y(:,kx,ky,kz);
%                     r2s(kx,ky,kz) = b(2);
%                     m0(kx,ky,kz) = exp(b(1));
%                 end
%             end
%         end
%     end
% else
%     for kz=1:nz
%         for ky=1:ny
%             for kx=1:nx
%                 if mask(kx,ky,kz)==1
%                     b = x\y(:,kx,ky,kz);
%                     r2s(kx,ky,kz) = b(2);
%                     m0(kx,ky,kz) = exp(b(1));
%                 end
%             end
%         end
%     end
% end
r2s = SetImgRange(r2s,ranger2s);
t2s = 1./r2s;
t2s = SetImgRange(t2s,ranget2s);
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
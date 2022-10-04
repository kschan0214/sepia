%% function [r2s, M0final] = R2starmapping_ARLO_complex(img,te)
%
% Description: Fast monoexponential fitting of R2* by auto-regression
% echo spacing should be the same across time series
% Assuming time series in last dimension
% M0 map is calculated by weighted sum method
% ref: Pei et al. MRM 73:843?850(2015)
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 11 October, 2016
% Date last modified: 
%
function [r2s, M0final] = R2starmapping_ARLO_complex(img,te)
% get dimension of data
matrixSize = size(img);
ndim = ndims(img);

% reshape image
img = reshape(img,[numel(img)/matrixSize(end),matrixSize(end)]);
% Remove the effect of time varying phase change and initial phase
img2 = zeros(numel(img)/matrixSize(end),matrixSize(end)-1);
for kt=1:length(te)-1
    img2(:,kt) = conj(img(:,kt)).*img(:,kt+1);
end
img2 = reshape(img2,[matrixSize(1:end-1),matrixSize(end)-1]);

dt = te(2)-te(1);
dx = 1:length(te)-1;

u = ARLO(img2,dx);

% u = 1/(2*R2s*dt)
r2s = 1./(real(u)*(2*dt));
r2s(r2s<0)=0;
r2s(r2s>5000)=5000;
r2s(isnan(r2s))=0;
r2s(isinf(r2s))=5000;

% compute M0map
% calculate weighting for each echo for combination
img = abs(reshape(img,[matrixSize(1:end-1),matrixSize(end)]));
total_M = sum(img,ndim);
w = bsxfun(@rdivide,img,total_M);
w(isnan(w))=0;
temat = permute(repmat(te(:),[1,matrixSize(1:end-1)]),[2:ndim,1]);
r2s_inv = bsxfun(@times,r2s,temat);
M0tmp = img.*r2s_inv;
% weighted cimbination based on signal intensity
M0final = sum(w.*M0tmp,ndim);

end
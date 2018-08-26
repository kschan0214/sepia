%% function m0 = ComputeM0GivenR2star(r2s,te,img,method)
%
% Input
% --------------
%   r2s         : R2* map
%   te          : echo times
%   img         : multiecho images
%   method      : '1stecho','weighted','average'
%
% Output
% --------------
%   m0          : T1-weighted images
%
% Description: compute T1-wighted images given R2* and TE
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 30 June 2017
% Date last modified:
%
%
function m0 = ComputeM0GivenR2star(r2s,te,img,method)
% get dimension of data
matrixSize = size(img);
ndim = ndims(img);

switch lower(method)
    case '1stecho'
        m0 = img(1:(numel(img)/matrixSize(end)))'.*exp(r2s(:)*te(1));
        if numel(m0) ~=1
            m0 = reshape(m0,matrixSize(1:end-1));
        end
    case 'weighted'
        % calculate weighting for each echo for combination
        total_M = sum(img,ndim);
        w = bsxfun(@rdivide,img,total_M);
        w(isnan(w))=0;
        teMat = permute(repmat(te(:),[1,matrixSize(1:end-1)]),[2:ndim,1]);
        r2sMat = bsxfun(@times,r2s,teMat);
        m0tmp = img.*exp(r2sMat);
        % weighted cimbination based on signal intensity
        m0 = sum(w.*m0tmp,ndim);
    case 'average'
        teMat = permute(repmat(te(:),[1,matrixSize(1:end-1)]),[2:ndim,1]);
        r2sMat = bsxfun(@times,r2s,teMat);
        m0tmp = img.*exp(r2sMat);
        % weighted cimbination based on signal intensity
        m0 = mean(m0tmp,ndim);
end
m0(isnan(m0)) = 0;
m0(isinf(m0)) = 0;
end
% [totalField, N_std, headerAndExtraData] = Wrapper_EchoCombine_OptimumWeight(fieldMap,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)
%
% Input
% --------------
% fieldMap      : original single-/multi-echo wrapped phase image, in rad
% mask          : signal mask
% matrixSize    : size of the input image
% voxelSize     : spatial resolution of each dimension of the data, in mm
% algorParam    : structure contains fields with algorithm-specific parameter(s)
% headerAndExtraData : structure contains extra header info/data for the algorithm
%
% Output
% --------------
% totalField    : unwrapped total field, in radHz
% N_std         : noise standard deviation in the field map
% fieldmapUnwrapAllEchoes : unwrapped echo phase, in rad
%
% Description: Wrapper to perform phase echo combination
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 22 June 2021 (v1.0)
% Date modified:
%
%
function [totalField, N_std, headerAndExtraData] = Wrapper_EchoCombine_OptimumWeight(fieldMap,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)

% get some data from headerAndExtraData
magn	= double(headerAndExtraData.magn);
TE      = headerAndExtraData.te;

% initiate empty variable
fieldmapUnwrapAllEchoes = [];

% check if the input data multiecho or not
isMultiecho = size(fieldMap,4) > 1;
% find the centre of mass
% pos=round(centerofmass(magn(:,:,:,1)));

if isMultiecho
    
%%%%%%%%%%%%%%%%%%%%%%%% Multi-echo %%%%%%%%%%%%%%%%%%%%%%%%
dims = size(fieldMap);
dims(4) = dims(4) - 1;
    
%%%%%%%%%%%%%%%%%%%%%%%% Step 1: unwrap echo phase %%%%%%%%%%%%%%%%%%%%%%%%
    % compute wrapped phase shift between successive echoes & unwrap each echo phase shift
    phaseShiftUnwrapAllEchoes = zeros(dims, 'like',fieldMap);
    for k = 1:dims(4)
        fprintf('Unwrapping #%i echo shift...\n',k);
        tmp	= angle(exp(1i*fieldMap(:,:,:,k+1))./exp(1i*fieldMap(:,:,:,k)));
        tmp = UnwrapPhaseMacro(tmp,mask,matrixSize,voxelSize,algorParam,headerAndExtraData);
%         tmp2(:,:,:,k) = tmp-round(tmp(pos(1),pos(2),pos(3))/(2*pi))*2*pi;
        phaseShiftUnwrapAllEchoes(:,:,:,k) = tmp-round(mean(tmp( mask == 1))/(2*pi))*2*pi;
    end
    % get phase accumulation over all echoes
    phaseShiftUnwrapAllEchoes = cumsum(phaseShiftUnwrapAllEchoes,4);   
    
    % get the unwrapped phase accumulation across echoes
    % unwrap first echo
    tmp                     = UnwrapPhaseMacro(fieldMap(:,:,:,1),mask,matrixSize,voxelSize,algorParam,headerAndExtraData);
    tmp                     = tmp-round(mean(tmp( mask == 1))/(2*pi))*2*pi;
    fieldmapUnwrapAllEchoes = cat(4,tmp,phaseShiftUnwrapAllEchoes + tmp);
%     fieldmapUnwrapAllEchoes = cumsum(fieldmapUnwrapAllEchoes,4);
    
%%%%%%%%%%%%%%%%%%%%%%%% Step 2: Compute weights %%%%%%%%%%%%%%%%%%%%%%%%
    % Robinson et al. 2017 NMR Biomed Appendix A2
    N_std = zeros(dims, 'like',fieldMap);   % fieldmap SD
    for k=1:dims(4)
        N_std(:,:,:,k) = 1./(TE(k+1)-TE(1)) ...
            * sqrt((magn(:,:,:,1).^2+magn(:,:,:,k+1).^2)./((magn(:,:,:,1).*magn(:,:,:,k+1)).^2));
    end
    % weights are inverse of the field map variance
    weight = bsxfun(@rdivide,1./(N_std.^2),sum(1./(N_std.^2),4));
    
    % standard deviation of field map from weighted avearging
    N_std               = sqrt(sum(weight.^2 .* N_std.^2,4));    % sqrt(Weighted variance) = SD
    N_std(isnan(N_std))	= 0;
    N_std(isinf(N_std))	= 0;
    N_std               = N_std./norm(N_std(mask~=0));
%     totalFieldSD                        = totalFieldSD./norm(rmoutliers(totalFieldSD(mask~=0)));

%%%%%%%%%%%%%%%%%%%%%%%% Step 3: Weighted average %%%%%%%%%%%%%%%%%%%%%%%%
% 20210803: use for-loop to reduce memory
totalField = zeros(dims(1:3), 'like',fieldMap);
for k = 1:dims(4)
    % Weighted average of unwrapped phase shift
    totalField = totalField + weight(:,:,:,k).*(phaseShiftUnwrapAllEchoes(:,:,:,k)/(TE(k+1)-TE(1)));
end
totalField(isnan(totalField)) = 0;
totalField(isinf(totalField)) = 0;

else

%%%%%%%%%%%%%%%%%%%%%%%% Single echo %%%%%%%%%%%%%%%%%%%%%%%%
    totalField	= UnwrapPhaseMacro(fieldMap,mask,matrixSize,voxelSize,algorParam,headerAndExtraData);
    totalField 	= totalField-round(mean(totalField( mask == 1))/(2*pi))*2*pi;
    totalField	= totalField/(TE(1));
    N_std            	= 1./magn;
    N_std(isnan(N_std))	= 0;
    N_std(isinf(N_std))	= 0;
    N_std               = N_std./norm(N_std(mask~=0));
%     totalFieldSD                        = totalFieldSD./norm(rmoutliers(totalFieldSD(mask~=0)));
end


% apply mask
totalField              = bsxfun(@times,totalField,mask);
N_std                   = bsxfun(@times,N_std,mask);
fieldmapUnwrapAllEchoes = bsxfun(@times,fieldmapUnwrapAllEchoes,mask);

if ~isempty(fieldmapUnwrapAllEchoes)
    headerAndExtraData.fieldmapUnwrapAllEchoes = fieldmapUnwrapAllEchoes;
end

end

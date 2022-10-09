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
% Date modified: 3 Aug 2021 (v1.0)
% Date modified: 13 Aug 2021 (v1.0)
%
%
function [totalField, N_std, headerAndExtraData] = Wrapper_EchoCombine_OptimumWeight(fieldMap,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)

% get extra data such as magnitude/weights/B0 direction/TE/etc.
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
TE      = headerAndExtraData.sepia_header.TE;
magn    = get_variable_from_headerAndExtraData(headerAndExtraData, 'magnitude', matrixSize);


% % initiate empty variable
% fieldmapUnwrapAllEchoes = [];

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
    %         fieldMap(:,:,:,k+1) = tmp-round(mean(tmp( mask == 1))/(2*pi))*2*pi;
    end
    % get phase accumulation over all echoes
    phaseShiftUnwrapAllEchoes = cumsum(phaseShiftUnwrapAllEchoes,4);   
    %     fieldMap(:,:,:,2:end) = cumsum(fieldMap(:,:,:,2:end),4); 

    % get the unwrapped phase accumulation across echoes
    % unwrap first echo
    tmp                     = UnwrapPhaseMacro(fieldMap(:,:,:,1),mask,matrixSize,voxelSize,algorParam,headerAndExtraData);
    tmp                     = tmp-round(mean(tmp( mask == 1))/(2*pi))*2*pi;
    %     fieldMap = cat(4,tmp,fieldMap(:,:,:,2:end) + tmp);
    % fieldmapUnwrapAllEchoes = cat(4,tmp,phaseShiftUnwrapAllEchoes + tmp);
    fieldMap = cat(4,tmp,phaseShiftUnwrapAllEchoes + tmp);
    clear tmp phaseShiftUnwrapAllEchoes % release memory

    %%%%%%%%%%%%%%%%%%%%%%%% Step 2: Compute weights %%%%%%%%%%%%%%%%%%%%%%%%
    % Robinson et al. 2017 NMR Biomed Appendix A2
    [weight, N_std]= compute_optimum_weighting_combining_phase_difference(magn, TE);
    clear magn % release memory

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
%         totalField = totalField + weight(:,:,:,k).*(phaseShiftUnwrapAllEchoes(:,:,:,k)/(TE(k+1)-TE(1)));
        totalField = totalField + weight(:,:,:,k).*((fieldMap(:,:,:,k+1) - fieldMap(:,:,:,1))/(TE(k+1)-TE(1)));
    end
    totalField(isnan(totalField)) = 0;
    totalField(isinf(totalField)) = 0;
    
    fieldMap = bsxfun(@times,fieldMap,mask);
    headerAndExtraData.fieldmapUnwrapAllEchoes = fieldMap; % use headerAndExtraData to pass fieldmapUnwrapAllEchoes back

else

    %%%%%%%%%%%%%%%%%%%%%%%% Single echo %%%%%%%%%%%%%%%%%%%%%%%%
    N_std            	= 1./magn;
    N_std(isnan(N_std))	= 0;
    N_std(isinf(N_std))	= 0;
    N_std               = N_std./norm(N_std(mask~=0));
    clear magn
    
    totalField	= UnwrapPhaseMacro(fieldMap,mask,matrixSize,voxelSize,algorParam,headerAndExtraData);
    totalField 	= totalField-round(mean(totalField( mask == 1))/(2*pi))*2*pi;
    totalField	= totalField/(TE(1));
    %     totalFieldSD                        = totalFieldSD./norm(rmoutliers(totalFieldSD(mask~=0)));
end


% apply mask
totalField              = bsxfun(@times,totalField,mask);
N_std                   = bsxfun(@times,N_std,mask);

end

%% [totalField,N_std,fieldmapUnwrapAllEchoes] = estimateTotalField(fieldMap,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)
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
% totalField    : unwrapped total field, in Hz
% N_std         : noise standard deviation in the field map
% fieldmapUnwrapAllEchoes : unwrapped echo phase, in rad
%
% Description: This is a wrapper function to access individual echo combination for phase
%              algorithms for SEPIA (default: 'MEDI non-linear')
%
% This code is modified from the T2starAndFieldCalc.m from Jose P. Marques
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 31 May, 2017
% Date modified: 27 February 2018
% Date modified: 16 September 2018
% Date modified: 24 may 2019
% Date modified: 9 March 2020 (v0.8.0)
%
function [totalField,N_std,fieldmapUnwrapAllEchoes] = estimateTotalField(fieldMap,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)
sepia_universal_variables;

matrixSize  = double(matrixSize(:).');
voxelSize   = double(voxelSize(:).');
fieldmapUnwrapAllEchoes = [];

algorParam	= check_and_set_SEPIA_algorithm_default(algorParam);
echoCombine	= algorParam.unwrap.echoCombMethod;
unit        = algorParam.unwrap.unit;

headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
TE = headerAndExtraData.te;
dt = headerAndExtraData.delta_TE;
if ~isempty(headerAndExtraData.magn)
    magn = double(headerAndExtraData.magn);
else
    magn = repmat(mask,1,1,1,size(fieldMap,4));
end

%% ensure all variables have the same data type
fieldMap	= double(fieldMap);
mask       	= double(mask);

% find the centre of mass
pos=round(centerofmass(magn(:,:,:,1)));

disp('--------------------');
disp('Total field recovery');
disp('--------------------');
disp('Calculating field map...');

fprintf('Temporal phase unwrapping: %s \n',echoCombine);

% if numel(TE) < 3
%     echoCombine = methodEchoCombineName{2};
%     warning('Number of echoes is less than 3. Using Optimum Weight method to combine echo phase instead.');
% end

%% Core
switch echoCombine
    case methodEchoCombineName{1}
        if size(fieldMap,4) > 1
        %%%%%%%%%%%%%%%%%%%%%%%% Multi-echo %%%%%%%%%%%%%%%%%%%%%%%%
            % compute wrapped phase shift between successive echoes
            fieldMapEchoTemp = angle(exp(1i*fieldMap(:,:,:,2:end))./exp(1i*fieldMap(:,:,:,1:end-1)));
            
            % unwrap each echo phase shift
            tmp2             = zeros(size(fieldMapEchoTemp),'like',fieldMap);
            for k = 1:size(fieldMapEchoTemp,4)
                fprintf('Unwrapping #%i echo shift...\n',k);
                tmp           = UnwrapPhaseMacro(fieldMapEchoTemp(:,:,:,k),mask,matrixSize,voxelSize,algorParam,headerAndExtraData);
                tmp2(:,:,:,k) = tmp-round(tmp(pos(1),pos(2),pos(3))/(2*pi))*2*pi;
            end
            
            % get phase accumulation over all echoes
            phaseShiftUnwrapAllEchoes=cumsum(tmp2,4);    

            % compute unwrapped phase shift between successive echoes
            fieldMapUnwrap = zeros(size(phaseShiftUnwrapAllEchoes),'like',fieldMap);
            for k = 1:size(phaseShiftUnwrapAllEchoes,4)
                fieldMapUnwrap(:,:,:,k) = phaseShiftUnwrapAllEchoes(:,:,:,k)/(TE(k+1)-TE(1));
            end

            % Robinson et al. 2017 NMR Biomed Appendix A2
            fieldMapSD = zeros(size(phaseShiftUnwrapAllEchoes),'like',fieldMap);
            for k=1:size(phaseShiftUnwrapAllEchoes,4)
                fieldMapSD(:,:,:,k) = 1./(TE(k+1)-TE(1)) ...
                    * sqrt((magn(:,:,:,1).^2+magn(:,:,:,k+1).^2)./((magn(:,:,:,1).*magn(:,:,:,k+1)).^2));
            end
            
            % weights are inverse of the field map variance
            weight = bsxfun(@rdivide,1/(fieldMapSD.^2),sum(1/(fieldMapSD.^2),4));
            
            % Weighted average of unwrapped phase shift
            totalField                    = sum(fieldMapUnwrap .* weight,4);
            totalField(isnan(totalField)) = 0;
            totalField(isinf(totalField)) = 0;
            
            % standard deviation of field map from weighted avearging
            totalFieldVariance                  = sum(weight.^2 .* fieldMapSD.^2,4);
            totalFieldSD                        = sqrt(totalFieldVariance);
            totalFieldSD(isnan(totalFieldSD))   = 0;
            totalFieldSD(isinf(totalFieldSD))   = 0;
            totalFieldSD                        = totalFieldSD./norm(totalFieldSD(mask~=0));
            
            % get the unwrapped phase accumulation across echoes
            % unwrap first echo
            tmp                     = UnwrapPhaseMacro(fieldMap(:,:,:,1),mask,matrixSize,voxelSize,algorParam,headerAndExtraData);
            tmp                     = tmp-round(tmp(pos(1),pos(2),pos(3))/(2*pi))*2*pi;
            fieldmapUnwrapAllEchoes = cat(4,tmp,tmp2);
            fieldmapUnwrapAllEchoes = cumsum(fieldmapUnwrapAllEchoes,4);
            
        else
            
        %%%%%%%%%%%%%%%%%%%%%%%% Single echo %%%%%%%%%%%%%%%%%%%%%%%%
            tmp                                 = UnwrapPhaseMacro(fieldMap,mask,matrixSize,voxelSize,algorParam,headerAndExtraData);
            tmp                                 = tmp-round(tmp(pos(1),pos(2),pos(3))/(2*pi))*2*pi;
            totalField                          = tmp/(TE(1));
            totalFieldSD                        = 1./magn;
            totalFieldSD(isnan(totalFieldSD))   = 0;
            totalFieldSD(isinf(totalFieldSD))   = 0;
            totalFieldSD                        = totalFieldSD./norm(totalFieldSD(mask~=0));
        end
        
        N_std = totalFieldSD;
        
    case methodEchoCombineName{2}
        sepia_addpath('MEDI');
        if numel(TE)>3 && ((TE(2)-TE(1))-(TE(3)-TE(2))>1e-5)
            % Estimate the frequency offset in each of the voxel using a complex
            % fitting (uneven echo spacing)
            [iFreq_raw, N_std] = Fit_ppm_complex_TE(magn.*exp(-1i*fieldMap),TE);
        else
            % Estimate the frequency offset in each of the voxel using a complex
            % fitting (even echo spacing)
            [iFreq_raw, N_std] = Fit_ppm_complex(magn.*exp(-1i*fieldMap));
        end

        % Compute magnitude image
        headerAndExtraData.magn = sqrt(sum(abs(magn).^2,4));
    
        % Spatial phase unwrapping
        totalField = UnwrapPhaseMacro(iFreq_raw,mask,matrixSize,voxelSize,algorParam,headerAndExtraData);
                    
        % use the centre of mass as reference phase
        totalField = totalField-round(totalField(pos(1),pos(2),pos(3))/(2*pi))*2*pi;
        
        % convert rad to radHz
        totalField = totalField / dt;
        
    case methodEchoCombineName{3}
        sepia_addpath('MEDI');
        % Estimate the frequency offset in each of the voxel using a complex
        % fitting (even echo spacing)
        [iFreq_raw, N_std] = Fit_ppm_complex_bipolar(magn.*exp(-1i*fieldMap));

        % Compute magnitude image
        headerAndExtraData.magn = sqrt(sum(abs(magn).^2,4));

        % Spatial phase unwrapping
        totalField = UnwrapPhaseMacro(iFreq_raw,mask,matrixSize,voxelSize,algorParam,headerAndExtraData);
                    
        % use the centre of mass as reference phase
        totalField = totalField-round(totalField(pos(1),pos(2),pos(3))/(2*pi))*2*pi;
        
        % convert rad to radHz
        totalField = totalField / dt;
        
end

disp(['The resulting field map with the following unit: ' unit]);
switch lower(unit)
    case 'ppm'
        totalField = (totalField/(2*pi))/(fieldStrength*gyro);
%         tmp2 = (tmp2/(2*pi))/(fieldStrength*gamma);
    case 'rad'
        totalField = totalField*dt;
%         tmp2 = tmp2*dt;
    case 'hz'
        totalField = totalField/(2*pi);
%         tmp2 = tmp2/(2*pi);
    case 'radhz'
    otherwise
        disp(['Input unit is invalid. radHz is used instead.']);
end

N_std = real(N_std);

end

%% find the centre of mass
function coord=centerofmass(data)
data=abs(data);
dims=size(data);
    for k=1:length(dims)
    %     datatemp=permute(data,[k ]);
    dimsvect=ones([1, length(dims)]);
    dimsvect(k)=dims(k);
    temp=bsxfun(@times,(data),reshape(1:dims(k),dimsvect));
    coord(k)=sum(temp(:))./sum(data(:));
    end
end

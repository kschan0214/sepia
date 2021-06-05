%% [totalField, N_std, headerAndExtraData] = Wrapper_EchoCombine_MEDInonlinearfit(fieldMap,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)
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
% Date created: 5 June 2021 (v1.0)
% Date last modified:
%
%
function [totalField, N_std, headerAndExtraData] = Wrapper_EchoCombine_MEDInonlinearfit(fieldMap,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)


% get some data from headerAndExtraData
magn	= double(headerAndExtraData.magn);
TE      = headerAndExtraData.te;
dt      = headerAndExtraData.delta_TE;

% find the centre of mass
pos     = round(centerofmass(magn(:,:,:,1)));

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

% temporary usage
headerAndExtraData.magn = sqrt(sum(abs(magn).^2,4));

% Spatial phase unwrapping
totalField = UnwrapPhaseMacro(iFreq_raw,mask,matrixSize,voxelSize,algorParam,headerAndExtraData);

% use the centre of mass as reference phase
totalField = totalField-round(totalField(pos(1),pos(2),pos(3))/(2*pi))*2*pi;

% convert rad to radHz
totalField = totalField / dt;
        
% return the original data
headerAndExtraData.magn = magn;

end


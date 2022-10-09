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
% Date modified: 5 June 2021 (v1.0)
% Date modified: 12 September 2022 (v1.1)
%
function [totalField,N_std,fieldmapUnwrapAllEchoes,mask] = estimateTotalField(fieldMap,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)
sepia_universal_variables;

matrixSize  = double(matrixSize(:).');
voxelSize   = double(voxelSize(:).');
fieldmapUnwrapAllEchoes = [];

algorParam	= check_and_set_SEPIA_algorithm_default(algorParam);
echoCombine	= algorParam.unwrap.echoCombMethod;
unit        = algorParam.unwrap.unit;

headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
dt = headerAndExtraData.sepia_header.delta_TE;

if isempty(headerAndExtraData.availableFileList.magnitude) && isempty(headerAndExtraData.magnitude)
    headerAndExtraData.magn = repmat(mask,1,1,1,size(fieldMap,4));
end

%% ensure all variables have the same data type
fieldMap	= double(fieldMap);
mask       	= double(mask);

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
for k = 1:length(wrapper_EchoCombine_function)
    if strcmpi(echoCombine,methodEchoCombineName{k})
        [totalField, N_std, headerAndExtraData] = feval(wrapper_EchoCombine_function{k},fieldMap,mask,matrixSize,voxelSize,algorParam,headerAndExtraData);
    end
end

if isfield(headerAndExtraData,'fieldmapUnwrapAllEchoes')
    fieldmapUnwrapAllEchoes = headerAndExtraData.fieldmapUnwrapAllEchoes;
    headerAndExtraData = rmfield(headerAndExtraData,'fieldmapUnwrapAllEchoes');
end

% in case the echo combine algorithm returns a different mask, e.g., ROMEO
if isfield(headerAndExtraData, 'mask')
    mask = headerAndExtraData.mask;
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

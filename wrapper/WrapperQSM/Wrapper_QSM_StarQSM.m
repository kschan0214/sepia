%% [chi] = Wrapper_QSM_StarQSM(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
%
% Input
% --------------
% localField    : local field map (tissue fields), in Hz
% mask          : signal mask
% matrixSize    : size of the input image
% voxelSize     : spatial resolution of each dimension of the data, in mm
% algorParam    : structure contains fields with algorithm-specific parameter(s)
% headerAndExtraData : structure contains extra header info/data for the algorithm
%
% Output
% --------------
% chi           : magnetic susceptibility map, in ppm
%
% Description: This is a wrapper function to access Star-QSM for SEPIA
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 8 March 2020
% Date modified: 13 August 2021 (v1.0)
%
%
function [chi] = Wrapper_QSM_StarQSM(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
sepia_universal_variables;

% get algorithm parameters
algorParam = check_and_set_algorithm_default(algorParam);
method     = algorParam.qsm.method;
padSize    = algorParam.qsm.padsize;

% get extra data such as magnitude/weights/B0 direction/TE/etc.
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
b0dir = headerAndExtraData.sepia_header.B0_dir;
b0    = headerAndExtraData.sepia_header.B0;
te    = headerAndExtraData.sepia_header.delta_TE;

% add path
sepia_addpath('STISuite');

%% Display algorithm parameters
disp('The following parameter is being used...');
disp(['Pas size = [' num2str(padSize(1)) ', ' num2str(padSize(2)) ', ' num2str(padSize(3)) ']']);

%% main
% Unlike the iLSQR implementation, the order of local field map
% values will affect the Star-QSM result, i.e. 
% chi = method(localField,...) ~= method(localField*C,...)/C, where
% C is a constant. Lower order of local field magnitude will more 
% likely produce chi map with streaking artefact. 
% In the STI_Templates.m example, Star-QSM expecting local field in 
% the unit of rad. However, value of the field map in rad will 
% vary with echo time. Therefore, data acquired with short
% TE will be prone to streaking artefact. To mitigate this
% potential problem, local field map is converted from Hz to radHz
% here and the resulting chi will be normalised by the same factor 
% 
localField = localField * 2*pi;

chi = QSM_star(localField,mask,'TE',te*1e3,'B0',b0,'H',b0dir,'padsize',padSize,'voxelsize',voxelSize);

chi = chi * te;
        

end

%% set default parameter if not specified
function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;

try algorParam2.qsm.padSize = algorParam.qsm.padSize; catch; algorParam2.qsm.padSize = [12,12,12]; end

end
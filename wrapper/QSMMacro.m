%% function [chi] = QSMMacro(localField,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)
%
% Input
% --------------
% localField    : local field map (tissue field), in Hz
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
% Description: This is a wrapper function to access individual dipole field
%              inversion algorithms for SEPIA (default: 'TKD')
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 28 June 2017
% Date modified: 9 April 2018
% Date modified: 1 April 2019
% Date modified: 5 June 2019
% Date modified: 27 Feb 2020 (v0.8.0)
% Date modified: 16 August 2021 (v1.0)
% Date modified: 19 July 2025 (v1.3)
%
function [chi,mask_ref,chi_para,chi_dia] = QSMMacro(localField,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)

sepia_universal_variables;
methodQSMName = lower(methodQSMName);

voxelSize       = double(voxelSize(:).');
matrixSize      = double(matrixSize(:).');

algorParam          = check_and_set_SEPIA_algorithm_default(algorParam);
method              = algorParam.qsm.method;
reference_tissue    = algorParam.qsm.reference_tissue;

headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);

disp('---------------------------');
disp('Dipole field inversion step');
disp('---------------------------');

%% zero padding for odd number dimension
fprintf('Zero-padding data if the input images have odd number matrix size...');
% essential input
localField  = double(zeropad_odd_dimension(localField,'pre'));
mask        = double(zeropad_odd_dimension(mask,'pre'));
matrixSize_new = size(localField);

% additional input
if ~isempty(headerAndExtraData.weights)
    headerAndExtraData.weights = double(zeropad_odd_dimension(headerAndExtraData.weights,'pre'));
end
if ~isempty(headerAndExtraData.magnitude)
    headerAndExtraData.magnitude = double(zeropad_odd_dimension(headerAndExtraData.magnitude,'pre'));
end
if ~isempty(headerAndExtraData.initGuess)
    headerAndExtraData.initGuess = double(zeropad_odd_dimension(headerAndExtraData.initGuess,'pre'));
end   

fprintf('Done!\n');

%% reference tissue
switch reference_tissue
    case 'None'
        mask_ref = [];
        
    case 'Brain mask'
        mask_ref = mask;
        
    case 'CSF'
        if( isempty(headerAndExtraData.magnitude) && isempty(headerAndExtraData.availableFileList.magnitude))
            warning('Please specify a magnitude data (at least 3 echoes) if you want to use CSF as reference.');
            warning('No normalisation will be done on the susceptibility map in this instance.');
            mask_ref = [];
        else
            sepia_addpath('MEDI');
            magn        = get_variable_from_headerAndExtraData(headerAndExtraData, 'magnitude', matrixSize_new);
            if size(magn,4) < 3
                warning('Please specify a magnitude data (at least 3 echoes) if you want to use CSF as reference.');
                warning('No normalisation will be done on the susceptibility map in this instance.');
                mask_ref = [];
                clear magn
                
            else
                r2s         = arlo(headerAndExtraData.sepia_header.TE, magn);
                clear magn

                mask_ref    = extract_CSF(r2s,mask,voxelSize)>0;
                clear r2s
            end
        end
end

%% QSM algorithm
disp('Computing QSM map...');
disp(['The following QSM algorithm will be used: ' method]);

% General steps as follow in the wrapper function:
% 1. input unit converted for optimal performance (if neccessary)
% 2. main QSM algorithm
% 3. convert output unit to ppm
for k = 1:length(wrapper_QSM_function)
    if strcmpi(method,methodQSMName{k})
        nOut = nargout(wrapper_QSM_function{k});
        varargout = cell(1, max(nOut,1));
        [varargout{:}] = feval(wrapper_QSM_function{k},localField,mask,matrixSize_new,voxelSize,algorParam, headerAndExtraData);
        chi = varargout{1};
        if nOut > 1
            chi_para = varargout{2};
            chi_dia  = varargout{3};
        end
        % [chi,varargout] = feval(wrapper_QSM_function{k},localField,mask,matrixSize_new,voxelSize,algorParam, headerAndExtraData);
    end
end

% Post dipole inversion using HEIDI
if algorParam.qsm.isHEIDI
    [chi] = Wrapper_QSM_HEIDI4all(localField,chi,mask,matrixSize,voxelSize,algorParam, headerAndExtraData) ;
end

% remove zero padding 
chi = double(zeropad_odd_dimension(chi,'post',matrixSize));
if exist("chi_para",'var')
    chi_para = double(zeropad_odd_dimension(chi_para,'post',matrixSize));
else
    chi_para = [];
end
if exist("chi_dia",'var')
    chi_dia = double(zeropad_odd_dimension(chi_dia,'post',matrixSize));
else
    chi_dia = [];
end
if ~isempty(mask_ref)
    mask_ref = double(zeropad_odd_dimension(mask_ref,'post',matrixSize));
end
mask_update = chi ~= 0;

% referencing
if ~isempty(mask_ref)
    if strcmpi(method,'MEDI')
        if ~algorParam.qsm.isLambdaCSF          % not MEDI+0, MEDI+0 needs no referencing
            chi(mask_update) = chi(mask_update) - mean(chi(mask_ref>0)); 
        else
            warning('MEDI+0 already uses CSF as reference region in optimisation. No referencing is performed.');
            mask_ref = [];
        end
    else
            chi(mask_update) = chi(mask_update) - mean(chi(mask_ref>0));
            if ~isempty(chi_para)
                chi_para(mask_update) = chi_para(mask_update) - mean(chi_para(mask_ref>0));
            end
            if ~isempty(chi_dia)
                chi_dia(mask_update) = chi_dia(mask_update) - mean(chi_dia(mask_ref>0));
            end

    end
end

end


%% Wrapper_SWISMWI_SWI_2DHamming(headerAndExtraData,output_structure,algorParam)
%
% Input
% --------------
% headerAndExtraData : structure contains extra header info/data for the algorithm
% output_structure   : structure contains output file path and output NIfTI header
% algorParam         : structure contains fields with algorithm-specific parameter(s)
%
% Output
% --------------
%
% Description: This is a wrapper function to access FANSI for SEPIA
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 3 August 2022 (v1.1)
% Date modified: 
%
%
function Wrapper_SWISMWI_SWI_2DHanning(headerAndExtraData,output_structure,algorParam)
sepia_universal_variables;

% get algorithm parameters
algorParam = check_and_set_algorithm_default(algorParam);
filterSize  = algorParam.swismwi.filterSize;
thres       = algorParam.swismwi.threshold;
m           = algorParam.swismwi.m;
method      = algorParam.swismwi.method;
ismIP       = algorParam.swismwi.ismIP;
slice_mIP   = algorParam.swismwi.slice_mIP;
isPositive  = algorParam.swismwi.isPositive;
isNegative  = algorParam.swismwi.isNegative;

% get extra data such as magnitude/weights/B0 direction/TE/etc.
headerAndExtraData                      = check_and_set_SEPIA_header_data(headerAndExtraData);
% check if the phase image is in wrapped phase range
headerAndExtraData.availableFileList	= io_true_phase_value(headerAndExtraData.availableFileList, output_structure.outputFileFullPrefix);
% b0dir = headerAndExtraData.sepia_header.B0_dir;
% matrixSize    = headerAndExtraData.sepia_header.matrixSize;
phase = get_variable_from_headerAndExtraData(headerAndExtraData, 'phasechi');
magn  = get_variable_from_headerAndExtraData(headerAndExtraData, 'magnitude');

outputNiftiTemplate     = output_structure.outputNiftiTemplate;
outputFileFullPrefix    = output_structure.outputFileFullPrefix;
% add path
% sepia_addpath;

%% SWI
disp('Computing SWI based on 2D Hamming filtering...');

[pSWI,nSWI,swi_phase] = swi_2dhanning(magn,phase,filterSize,thres,m,method);

% minimum intensity projection
if ismIP
    pswi_mIP = zeros(size(swi_phase));
    nswi_mIP = zeros(size(swi_phase));
    for kz = 1:size(swi_phase,3) - slice_mIP
        pswi_mIP(:,:,kz,:) = min(pSWI(:,:,kz:kz+slice_mIP-1,:),[],3);
        nswi_mIP(:,:,kz,:) = min(nSWI(:,:,kz:kz+slice_mIP-1,:),[],3);
    end
end

%% save nii
disp('Saving SWI results...');
save_nii_quick(outputNiftiTemplate, swi_phase, [outputFileFullPrefix 'swi-phase.nii.gz']);
if isPositive
    save_nii_quick(outputNiftiTemplate, pSWI, [outputFileFullPrefix 'swi-positive.nii.gz']);
    if ismIP
        save_nii_quick(outputNiftiTemplate, pswi_mIP, [outputFileFullPrefix 'minIP-positive.nii.gz']);
    end
end
if isNegative
    save_nii_quick(outputNiftiTemplate,nSWI, [outputFileFullPrefix 'swi-negative.nii.gz']);
    if ismIP
        save_nii_quick(outputNiftiTemplate,nswi_mIP, [outputFileFullPrefix 'minIP-negative.nii.gz']);
    end
end

disp('Done!');

end

%% set default parameter if not specified
function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;

try algorParam2.swismwi.filterSize	= algorParam.swismwi.filterSize;	catch; algorParam2.swismwi.filterSize	= 12;	end
try algorParam2.swismwi.threshold   = algorParam.swismwi.threshold;     catch; algorParam2.swismwi.threshold    = pi;	end
try algorParam2.swismwi.m           = algorParam.swismwi.m;             catch; algorParam2.swismwi.m            = 4;	end
try algorParam2.swismwi.ismIP   	= algorParam.swismwi.ismIP;         catch; algorParam2.swismwi.ismIP        = true;	end
try algorParam2.swismwi.slice_mIP  	= algorParam.swismwi.slice_mIP; 	catch; algorParam2.swismwi.slice_mIP	= 4;	end
try algorParam2.swismwi.isPositive	= algorParam.swismwi.isPositive; 	catch; algorParam2.swismwi.isPositive	= true;	end
try algorParam2.swismwi.isNegative	= algorParam.swismwi.isNegative; 	catch; algorParam2.swismwi.isNegative	= false;end
try algorParam2.swismwi.method      = algorParam.swismwi.method;        catch; algorParam2.swismwi.method       = 'default';end

end

%% I/O: convert phase to radian unit if required
function availableFileList          = io_true_phase_value(availableFileList, outputFileFullPrefix)

% load phase image to check if the 
phaseNIFTI = load_untouch_nii(availableFileList.phasechi);
phaseIMG   = load_nii_img_only(availableFileList.phasechi);

if abs(max(phaseIMG(:))-pi)>0.1 || abs(min(phaseIMG(:))-(-pi))>0.1 % allow small differences possibly due to data stype conversion or DICOM digitisation
% if abs(max(phaseNIFTI.img(:))-pi)>0.1 || abs(min(phaseNIFTI.img(:))-(-pi))>0.1 % allow small differences possibly due to data stype conversion or DICOM digitisation
    
    outputPhaseRadian_filename = [outputFileFullPrefix 'part-phase_rad.nii.gz'];

    disp('Values of input phase map exceed the range of [-pi,pi]. DICOM value is assumed.')
    fprintf('Rescaling phase data from DICOM image value to wrapped radian unit...')
    phase = DICOM2Phase(phaseNIFTI);
    fprintf('Done.\n')

    fprintf('Saving phase images in unit of radian...');
    save_nii_quick(phaseNIFTI, phase, outputPhaseRadian_filename);
    fprintf('Done.\n')
    
    % update the phase data for subsequent processing
    availableFileList.phasechi = outputPhaseRadian_filename;
    
end

end
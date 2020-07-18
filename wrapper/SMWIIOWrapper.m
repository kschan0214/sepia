%% SMWIIOWrapper(input,output,algorParam)
%
% Input
% --------------
% input         :   input structure containing filenames  
% output        :   output directory that stores the output 
% algorParam    :   structure contains method and method specific parameters
%
% Output
% --------------
%
% Description: Compute susceptibility map-weighted images
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 14 April 2019
% Date modified:
%
%
function SMWIIOWrapper(input,output,algorParam)

%% define variables
prefix = 'sepia_'; 

%% Check if output directory exists 
output_index = strfind(output, filesep);
outputDir = output(1:output_index(end));
% get prefix
if ~isempty(output(output_index(end)+1:end))
    prefix = [output(output_index(end)+1:end) '_'];
end
% if the output directory does not exist then create the directory
if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end

%% Check and set default algorithm parameters
algorParam      = CheckAndSetDefault(algorParam);
thres           = algorParam.smwi.threshold;
m               = algorParam.smwi.m;
ismIP           = algorParam.smwi.ismIP;
slice_mIP       = algorParam.smwi.slice_mIP;
isParamagnetic  = algorParam.smwi.isParamagnetic;
isDiamagnetic   = algorParam.smwi.isDiamagnetic;
isCSFRef        = algorParam.smwi.isCSFRef;

%% Read input
disp('Reading data...');

inputNiftiList = input;
    
if ~isempty(inputNiftiList(1).name)
    inputPhaseNifti = load_untouch_nii([inputNiftiList(1).name]);
    % load true value from nifti
    qsm = load_nii_img_only([inputNiftiList(1).name]);
%     qsm = double(inputPhaseNifti.img);
    
    if size(qsm,4) >1
        error('QSM map should be 3D. Please check your input.');
    end
    
    disp('QSM map is loaded.')
else
    error('Please specify a QSM map.');
end

if ~isempty(inputNiftiList(2).name)
    inputMagnNifti = load_untouch_nii([inputNiftiList(2).name]);
    % load tru value from nifti
    magn = load_nii_img_only([inputNiftiList(2).name]);
%     magn = double(inputMagnNifti.img);
    disp('Magnitude data is loaded.');
else
    error('Please specify a magnitude data.');
end

outputNiftiTemplate = inputMagnNifti;

%% SWI
disp('Computing SMWI...');

[pSMWI, dSMWI] = smwi(magn,qsm,thres,m);

% minimum intensity projection
if ismIP
    psmwi_mIP = zeros(size(pSMWI));
    dsmwi_mIP = zeros(size(pSMWI));
    for kz = 1:size(pSMWI,3) - slice_mIP
        psmwi_mIP(:,:,kz,:) = min(pSMWI(:,:,kz:kz+slice_mIP-1,:),[],3);
        dsmwi_mIP(:,:,kz,:) = min(dSMWI(:,:,kz:kz+slice_mIP-1,:),[],3);
    end
end

%% save nii
disp('Saving SMWI results...');
if isParamagnetic
    save_nii_quick(outputNiftiTemplate,pSMWI, [outputDir filesep prefix 'smwi-paramagnetic.nii.gz']);
    if ismIP
        save_nii_quick(outputNiftiTemplate,psmwi_mIP, [outputDir filesep prefix 'smwi-mIP-paramagnetic.nii.gz']);
    end
end
if isDiamagnetic
    save_nii_quick(outputNiftiTemplate,dSMWI, [outputDir filesep prefix 'smwi-diamagnetic.nii.gz']);
    if ismIP
        save_nii_quick(outputNiftiTemplate,dsmwi_mIP, [outputDir filesep prefix 'smwi-mIP-diamagnetic.nii.gz']);
    end
end

disp('Done!');

end

function algorParam2 = CheckAndSetDefault(algorParam)
algorParam2 = algorParam;

try algorParam2.smwi.threshold      = algorParam.smwi.threshold;        catch; algorParam2.smwi.threshold       = 1;	end
try algorParam2.smwi.m            	= algorParam.smwi.m;                catch; algorParam2.smwi.m             	= 4;	end
try algorParam2.smwi.ismIP         	= algorParam.smwi.ismIP;            catch; algorParam2.smwi.ismIP          	= true; end
try algorParam2.smwi.slice_mIP  	= algorParam.smwi.slice_mIP;        catch; algorParam2.smwi.slice_mIP      	= 4;	end
try algorParam2.smwi.isParamagnetic	= algorParam.smwi.isParamagnetic;   catch; algorParam2.smwi.isParamagnetic	= true;	end
try algorParam2.smwi.isDiamagnetic	= algorParam.smwi.isDiamagnetic;    catch; algorParam2.smwi.isDiamagnetic	= false;end
try algorParam2.smwi.isCSFRef       = algorParam.smwi.isCSFRef;         catch; algorParam2.smwi.isCSFRef        = false;end

end
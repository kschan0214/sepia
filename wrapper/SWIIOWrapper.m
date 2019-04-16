%% SWIIOWrapper(input,output,algorParam)
%
% Input
% --------------
%
% Output
% --------------
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 14 April 2019
% Date last modified:
%
%
function SWIIOWrapper(input,output,algorParam)

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
algorParam  = CheckAndSetDefault(algorParam);
filterSize  = algorParam.swi.filterSize;
thres       = algorParam.swi.threshold;
m           = algorParam.swi.m;
method      = algorParam.swi.method;
ismIP       = algorParam.swi.ismIP;
slice_mIP   = algorParam.swi.slice_mIP;
isPositive  = algorParam.swi.isPositive;
isNegative  = algorParam.swi.isNegative;

%% Read input
disp('Reading data...');

inputNiftiList = input;
    
if ~isempty(inputNiftiList(1).name)
    inputPhaseNifti = load_nii_4sepia([inputNiftiList(1).name]);
    phase = double(inputPhaseNifti.img);
    % check whether phase data contains DICOM values or wrapped
    % phase value
    if max(phase(:))>1000
        disp('Converting phase data from DICOM image value to radian unit...')
        phase = DICOM2Phase(inputPhaseNifti);

        disp('Saving phase images in unit of radian...');
        save_nii_quick(inputPhaseNifti,phase, [outputDir filesep prefix 'phase-rad.nii.gz']);

    end
    disp('Phase data is loaded.')
else
    error('Please specify a phase data.');
end

if ~isempty(inputNiftiList(2).name)
    inputMagnNifti = load_nii_4sepia([inputNiftiList(2).name]);
    magn = double(inputMagnNifti.img);
    disp('Magnitude data is loaded.');
else
    error('Please specify a magnitude data.');
end

outputNiftiTemplate = inputMagnNifti;

%% SWI
disp('Computing SWI...');

[pSWI,nSWI,swi_phase] = swi(magn,phase,filterSize,thres,m,method);

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
save_nii_quick(outputNiftiTemplate,swi_phase, [outputDir filesep prefix 'swi-phase.nii.gz']);
if isPositive
    save_nii_quick(outputNiftiTemplate,pSWI, [outputDir filesep prefix 'swi-positive.nii.gz']);
    if ismIP
        save_nii_quick(outputNiftiTemplate,pswi_mIP, [outputDir filesep prefix 'swi-mIP-positive.nii.gz']);
    end
end
if isNegative
    save_nii_quick(outputNiftiTemplate,nSWI, [outputDir filesep prefix 'swi-negative.nii.gz']);
    if ismIP
        save_nii_quick(outputNiftiTemplate,nswi_mIP, [outputDir filesep prefix 'swi-mIP-negative.nii.gz']);
    end
end

disp('Done!');

end

function algorParam2 = CheckAndSetDefault(algorParam)
algorParam2 = algorParam;

try algorParam2.swi.filterSize	= algorParam.swi.filterSize;	catch; algorParam2.swi.filterSize	= 12;	end
try algorParam2.swi.threshold   = algorParam.swi.threshold;     catch; algorParam2.swi.threshold    = pi;	end
try algorParam2.swi.m           = algorParam.swi.m;             catch; algorParam2.swi.m            = 4;	end
try algorParam2.swi.ismIP   	= algorParam.swi.ismIP;         catch; algorParam2.swi.ismIP        = true;	end
try algorParam2.swi.slice_mIP  	= algorParam.swi.slice_mIP; 	catch; algorParam2.swi.slice_mIP	= 4;	end
try algorParam2.swi.isPositive	= algorParam.swi.isPositive; 	catch; algorParam2.swi.isPositive	= true;	end
try algorParam2.swi.isNegative	= algorParam.swi.isNegative; 	catch; algorParam2.swi.isNegative	= false;end
try algorParam2.swi.method      = algorParam.swi.method;        catch; algorParam2.swi.method       = 'default';end

end
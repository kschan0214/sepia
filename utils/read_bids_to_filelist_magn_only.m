%% inputNIFTIList = read_bids_to_filelist(inputDir,outputPrefix)
%
% Input
% --------------
% inputDir      : input directory with BIDS
% outputPrefix  : full output basename
%
% Output
% --------------
% inputNIFTIList: file list for SEPIA
%
% Description: This script  (1) reads BIDS compatible input directory to a file list that can be used in SEPIA
%                           (2) converts separated echo images into 4D file (if needed)
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 11 August 2021 (v1.0)
% Date modified:
%
%
function inputNIFTIList = read_bids_to_filelist_magn_only(inputDir,outputPrefix)

disp('############################################');
disp('# Checking input directory for BIDS format #');
disp('############################################');

% default intput nifti list
inputNIFTIList(1).name = [outputPrefix 'part-phase.nii.gz'];
inputNIFTIList(2).name = [outputPrefix 'part-mag.nii.gz'];
inputNIFTIList(3).name = [];
inputNIFTIList(4).name = [outputPrefix 'header.mat'];

% check for NIFTI files with magnitude label
pattern = 'part-mag';
ext     = 'nii';
[magFile,magNumFiles] = get_filename_in_directory(inputDir,pattern,ext);
% check for JSON files with magnitude label
pattern = 'part-mag';
ext     = 'json';
[jsonFile,jsonNumFiles] = get_filename_in_directory(inputDir,pattern,ext);

% error if one of the required files cannot be found
if magNumFiles == 0 
    error('No magnitude file is found. For BIDS compatibility, make sure the magnitude NIFTI has the key ''part-mag''.');
end
if jsonNumFiles == 0 
    error('No JSON file is found. For BIDS compatibility, make sure the JSON has the key ''part-mag''.');
end

% preprocess files
if magNumFiles == 1 && jsonNumFiles == 1
    
    % Single-echo or single 4D volume route
    fprintf('One magnitude image is found:  %s \n', magFile(1).name);

    % magnitude
    inputNIFTIList(2).name = magFile(1).name;
    
else
    
    % Multi-echo, multiple volumes route
    % Step 1
    % 1.1 magnitude files
    disp('Multiple magnitude files are detected.');
    % get files that have 'echo' key
    magFile = validate_multiecho_key(magFile);
    % check if the files have consistent filename with BIDS specification
    validate_multiecho_single_acquisition(magFile);
    % update number of files
    magNumFiles = length(magFile);

    % 1.3 JSON files
    disp('Multiple JSON files are detected.');
    % get files that have 'echo' key
    jsonFile = validate_multiecho_key(jsonFile);
    % check if the files have consistent filename with BIDS specification
    validate_multiecho_single_acquisition(jsonFile);
    % update number of files
    jsonNumFiles = length(jsonFile);


    % magnitude 
    fprintf('Saving multi-echo magnitude data into a single volume...')
    isPhase = false;
    save_nifti_as_4d(magFile, inputNIFTIList(2).name, isPhase);
    fprintf('Done.\n')
    
end

% SEPIA header
save_sepia_header_from_bids(inputNIFTIList(2).name, jsonFile, outputPrefix);

end

%% Validate if the input filenames contains BIDS 'echo' key
function fileList = validate_multiecho_key(fileList)

fprintf('Checking if the files are multi-echo compatible...');

numFile = length(fileList);

% magnitude files
for k = numFile:-1:1
    if ~ContainName(fileList(k).name,'echo-')
        fileList(k) = [];
    end
end
% in case no file left then return error message
if isempty(fileList)
    disp('Failed!')
    error('No files contain BIDS key ''echo-'' in the filename.');
end

disp('Passed!');

end

%% Validate if the input filenames coming from a single acquisition
function validate_multiecho_single_acquisition(fileList)

fprintf('Checking if the files coming from a single acquisition...');

bids_basename = remove_bids_key(fileList(1).name, 'echo');

for k = 2:length(fileList)
    if ~strcmpi(bids_basename, remove_bids_key(fileList(k).name, 'echo'))
        disp('Failed!')
        error('It seems your BIDS directory contains multiple GRE acquisitions.\nFor example, \n%s and \n%s do not have consistent filenames other than the ''echo'' key.\n',fileList(1).name, fileList(k).name)
    end
end

disp('Passed!');

end

%% save multi-echo data from mutiple volume to a single image
function save_nifti_as_4d(fileList, outputFilename, isPhase)

if nargin < 3
    isPhase = false;
end

NumFiles = length(fileList);

numFileLoaded = 0;
img = [];
while numFileLoaded ~= NumFiles
    
    for k = 1:NumFiles
        
        if ContainName(fileList(k).name, ['echo-' num2str(numFileLoaded+1) '_']) || ContainName(fileList(k).name, ['echo-' num2str(numFileLoaded+1) '.'])
            nii = load_untouch_nii(fileList(k).name);
            
            if (abs(max(nii.img(:))-pi)>1e-4 || abs(min(nii.img(:))-(-pi))>1e-4) && isPhase % allow small differences possibly due to data stype conversion
                img	= cat(4, img, DICOM2Phase(nii));
            else
%                 img = cat(4, img, nii.img);
                % 20230307: bug fix applying rescale slope and intercept
                img = cat(4, img, load_nii_img_only(fileList(k).name));
            end
            
            numFileLoaded = numFileLoaded + 1;
        end
    end

end
save_nii_quick(nii,img,outputFilename);

end

%% save sepia header
function save_sepia_header_from_bids(niiFilename, jsonList, outpitPrefix)

input.nifti         = niiFilename;

for k = 1:length(jsonList)
    TEFileList{k} = jsonList(k).name;
end
input.TEFileList    = TEFileList;

save_sepia_header(input,[],outpitPrefix);

end


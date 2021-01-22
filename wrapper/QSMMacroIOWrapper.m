%% chi = QSMMacroIOWrapper(input,output,maskFullName,algorParam)
%
% Input
% --------------
% input         :   input directory contains NIfTI (*localfield*, *magn* and
%                   *fieldmapsd*) files or structure containing filenames  
% output        :   output directory that stores the output (susceptibility map)
% maskFullName  :   mask filename
% algorParam    :   structure contains method and method specific parameters
%
% Output
% --------------
% chi           : magnetic susceptibility map (in ppm)
%
% Description: This is a wrapper of QSMMacro.m which for NIfTI input/output
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 17 April 2018
% Date modified: 26 August 2018
% Date modified: 29 March 2019
% Date modified: 5 June 2019
% Date modified: 8 March 2020 (v0.8.0)
% Date modified: 21 Jan 2020 (v0.8.1)
%
%
function chi = QSMMacroIOWrapper(input,output,maskFullName,algorParam)
%% add general Path
sepia_addpath

sepia_universal_variables;

%% define variables
prefix = 'sepia_';
% make sure the input only load once (first one)
isLocalFieldLoad    = false;
isWeightLoad        = false;
isMagnLoad          = false; 

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

% display output info
fprintf('Output directory       : %s\n',outputDir);
fprintf('Output filename prefix : %s\n',prefix);

%% Setting up Input
disp('---------');
disp('Load data');
disp('---------');

% Step 1: check input for nifti files first
if isstruct(input)
    
    % Option 1: input are files
    inputNiftiList = input;
    
    % take the phase data directory as reference input directory 
    [inputDir,~,~] = fileparts(inputNiftiList(1).name);
    
else
    
    % Option 2: input is a directory
    disp('Searching input directory...');
    inputDir = input;
    
    % check and get filenames
    inputNiftiList  = struct();
    filePattern     = {'local-field','mag','weights','header'}; % don't change the order
    for  k = 1:length(filePattern)
        
        % get filename
        if k ~= 4   % NIFTI image input
            [file,numFiles] = get_filename_in_directory(inputDir,filePattern{k},'.nii');
        else        % SEPIA header
            [file,numFiles] = get_filename_in_directory(inputDir,filePattern{k},'.mat');
        end
        
        % actions given the number of files detected
        if numFiles == 1        % only one file -> get the name
            
            fprintf('One ''%s'' file is found: %s\n',filePattern{k},file.name);
            inputNiftiList(k).name = file.name;
            
        elseif numFiles == 0     % no file -> fatal error
            
            if k ~= 2 && k ~= 3 % essential file 'local-field' and header, i.e. k=1&4
                error(['No file with name containing string ''' filePattern{k} ''' is detected.']);
            else
                disp(['No file with name containing string ''' filePattern{k} ''' is detected.']);
            end
            
        else % multiple files -> fatal error
            error(['Multiple files with name containing string ''' filePattern{k} ''' are detected. Make sure the input directory should contain only one file with string ''' filePattern{k} '''.']);
        end
    end
end

%%%%%% Step 2: load data
% 2.1 Local field map  
if ~isempty(inputNiftiList(1).name)
    
    % load nifti structure for output template
    inputLocalFieldNifti = load_untouch_nii([inputNiftiList(1).name]);
    % load true value from NIfTI
    localField = load_nii_img_only([inputNiftiList(1).name]);
    isLocalFieldLoad = true;
    disp('Local field map is loaded.')
    
else
    error('Please specify a 3D local field map.');
end

% 2.2 magnitude data 
if ~isempty(inputNiftiList(2).name)
    
    % load true value from NIfTI
    magn = load_nii_img_only([inputNiftiList(2).name]);
    isMagnLoad = true;
    disp('Magnitude data is loaded.');
    
else
    disp('No magnitude data is loaded.');
    magn = ones(size(localField));
end

% 2.3 weights map 
if ~isempty(inputNiftiList(3).name)
    
    % load true value from NIfTI
    weights = load_nii_img_only([inputNiftiList(3).name]);
    isWeightLoad = true;
    disp('Weights data is loaded');
    
else
    disp('Default weighting method will be used for QSM.');
end
% 2.4 SEPIA header 
if ~isempty(inputNiftiList(4).name)
    load([inputNiftiList(4).name]);
    disp('Header data is loaded.');
else
    error('Please specify a header required by SEPIA.');
end

% if no magnitude or weighting map is loaded and MEDI is chosen -> fatal error
if ~isMagnLoad && strcmpi(algorParam.qsm.method,'MEDI')
    error('MEDI requires magnitude data. Please put the magnitude multi-echo data to input directory or use other algorithm');
end
if ~isWeightLoad && strcmpi(algorParam.qsm.method,'MEDI')
    error('MEDI requires a weighting map. Please put a (SNR) weighting map to input directory or use other algorithm');
end

% store the header the NIfTI files, all following results will have
% the same header
outputNiftiTemplate = inputLocalFieldNifti;


%%%%%% Step 3: validate input
% 3.1 Validate header information
validate_sepia_header;

% 3.2 Validate NIfTI input
disp('Validating input NIfTI files...')
if ~isequal(size(localField),matrixSize)
    erro('Input NIfTI files and SEPIA header do not have the same matrix size. Please check these files and/or remove the ''matrixSize'' variable from the SEPIA header.')
end
% check dimension of weights
if exist('weights','var')
    if ~isequal(size(weights),matrixSize)
        erro('Input NIfTI files and SEPIA header do not have the same matrix size. Please check these files and/or remove the ''matrixSize'' variable from the SEPIA header.')
    end
    if size(weights,4) > 1
        error('Input weighting map is 4D. SEPIA accepts weighting map to be 3D only.');
    end
end
disp('Input NIfTI files are valid.')

% display some header info
disp('----------------------');
disp('Basic data information');
disp('----------------------');
fprintf('Voxel size(x,y,z)   = %s mm x %s mm x %s mm\n' ,num2str(voxelSize(1)),num2str(voxelSize(2)),num2str(voxelSize(3)));
fprintf('Matrix size(x,y,z)  = %s x %s x %s\n'          ,num2str(matrixSize(1)),num2str(matrixSize(2)),num2str(matrixSize(3)));
fprintf('B0 direction(x,y,z) = [%s; %s; %s]\n'          ,num2str(B0_dir(1)),num2str(B0_dir(2)),num2str(B0_dir(3)));
fprintf('Field strength      = %s T\n'                  ,num2str(B0));
fprintf('Number of echoes    = %s\n'                    ,num2str(length(TE)));
fprintf('TE1/dTE             = %s/%s ms\n'              ,num2str(TE(1)*1e3),num2str(delta_TE*1e3));

% ensure variables are double
localField	= double(localField);
voxelSize   = double(voxelSize);
matrixSize  = double(matrixSize);
TE          = double(TE);
if exist('weights','var');  weights  = double(weights);	end
if exist('magn','var');     magn     = double(magn);   	end

%%%%%% Step 5: store some data to headerAndExtraData
% header
headerAndExtraData.b0           = B0;
headerAndExtraData.b0dir        = B0_dir;
headerAndExtraData.te           = TE;
headerAndExtraData.delta_TE     = delta_TE;
headerAndExtraData.CF           = CF;
headerAndExtraData.voxelSize    = voxelSize;
headerAndExtraData.matrixSize   = matrixSize;

if exist('magn','var');     headerAndExtraData.magn     = magn; end
if exist('weights','var');  headerAndExtraData.weights  = weights; end

clearvars inputLocalFieldNifti

%% get brain mask
disp('-----------');
disp('Signal mask');
disp('-----------');

% check if there is any mask specifically for dipole inversion first
maskList = dir(fullfile(inputDir, '*mask-qsm*'));
% if no final mask then just look for normal mask
if isempty(maskList)
    maskList = dir(fullfile(inputDir, '*mask*'));
end

% Scenario: No specified mask file + there is a file called mask in the input directory
if isempty(maskFullName) && ~isempty(maskList)
    
    fprintf('No mask file is specified but a mask file is found in the input directory: %s\n',fullfile(inputDir, maskList(1).name));
    disp('Trying to load the file as signal mask');
    
    maskFullName = fullfile(inputDir, maskList(1).name);
end

% Scenario: User provided a mask file or above scenario was satified
if ~isempty(maskFullName)
    
    % load mask file
    mask = load_nii_img_only(maskFullName) > 0;
    
    % make sure the mask has the same dimension as other input data
    if ~isequal(size(mask),matrixSize)
        disp('The file does not have the same dimension as other images.')
        % display error message 
        error('No mask file is loaded. Please specify a valid mask file or put it in the input directory.');
    else
        disp('Mask file is loaded.');
    end
    
else
    % display error message if nothing is found
    error('No mask file is found. Please specify your mask file or put it in the input directory.');
end

mask = double(mask);

if exist('weights','var');  headerAndExtraData.weights  = weights .* mask; end

% if ~isempty(maskFullName)
%     % Option 1: mask file is provided
%     mask = load_nii_img_only(maskFullName) > 0;
%     
% elseif ~isempty(maskList) 
%     % Option 2: input directory contains NIfTI file with name '*mask*'
%     fprintf('A mask file is found in the input directory: %s\n',fullfile(inputDir, maskList(1).name));
%     disp('Trying to load the file as signal mask');
%     mask = load_nii_img_only(fullfile(inputDir, maskList(1).name)) > 0;
%     
%     % make sure the mask has the same dimension as other input data
%     if ~isequal(size(mask),matrixSize)
%         disp('The file does not have the same dimension as other images.')
%         % display error message if nothing is found
%         error('No mask file is loaded. Pleasee specific your mask file or put it in the input directory.');
%     end
%     
%     disp('Mask file is loaded.');
%     
% else
%     % display error message if nothing is found
%     error('No mask file is loaded. Pleasee specific your mask file or put it in the input directory.');
%     
% end

%% Dipole inversion
% core of QSM
[chi,mask_ref] = QSMMacro(localField,mask,matrixSize,voxelSize,algorParam,headerAndExtraData);

% save results
fprintf('Saving susceptibility map...');

save_nii_quick(outputNiftiTemplate, chi, [outputDir filesep prefix 'QSM.nii.gz']);

if ~isempty(mask_ref)
    save_nii_quick(outputNiftiTemplate, mask_ref, [outputDir filesep prefix 'mask_reference_region.nii.gz']);
end

fprintf('Done!\n');

disp('Processing pipeline is completed!');

end

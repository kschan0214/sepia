%% save_sepia_header(input,userDefine,output)
%
% Usage:
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
% Date created: 23 May 2019
% Date last modified:
%
%
function save_sepia_header(input,userDefine,output)
% constant
gyro = 42.57747892;
prefix = 'sepia_';

% initiate all required header variables
B0_dir      = [];
B0          = [];
voxelSize   = [];
matrixSize  = [];
TE          = [];
delta_TE    = [];
CF          = [];

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

%% extract header info from input
if isstruct(input)
    %% Option 1: input are files
    % extract information from NIfTI
    
    % B0 direction
    inputNifti = load_untouch_nii(input.nifti);
    a = sqrt(1 - inputNifti.hdr.hist.quatern_b^2 - inputNifti.hdr.hist.quatern_c^2 - inputNifti.hdr.hist.quatern_d^2);
    rotmat_nifti=qGetR([a, inputNifti.hdr.hist.quatern_b,inputNifti.hdr.hist.quatern_c,inputNifti.hdr.hist.quatern_d]);
    B0_dir = rotmat_nifti \ [0;0;1];
    
    % voxel size
    voxelSize = inputNifti.hdr.dime.pixdim(2:4);
    
    % matrix size
    matrixSize = inputNifti.hdr.dime.dim(2:4);
    
    % TE
    if ~isempty(input.TEFileList)
        TE = readEchoTime(input.TEFileList);
    end
    
else
    % find all nifti files
    inputNiftiList = dir([input '/*.nii*']);
    
    if ~isempty(inputNiftiList)
        %% Option 2: input is a directory containing nifti and TE files
        % B0 direction
        inputNifti = load_untouch_nii(fullfile(input,inputNiftiList(1).name));
        a = sqrt(1 - inputNifti.hdr.hist.quatern_b^2 - inputNifti.hdr.hist.quatern_c^2 - inputNifti.hdr.hist.quatern_d^2);
        rotmat_nifti=qGetR([a, inputNifti.hdr.hist.quatern_b,inputNifti.hdr.hist.quatern_c,inputNifti.hdr.hist.quatern_d]);
        B0_dir = rotmat_nifti \ [0;0;1];

        % voxel size
        voxelSize = inputNifti.hdr.dime.pixdim(2:4);

        % matrix size
        matrixSize = inputNifti.hdr.dime.dim(2:4);

        inputTEList = dir([input '/*.txt*']);
        if ~isempty(inputTEList)
            TE = readEchoTime(inputTEList(1));
        else
            inputTEList = dir([input '/*.json*']);
            if ~isempty(inputTEList)
                
                for klist = 1:length(inputTEList)
                    inputTEListCell{klist} = fullfile(input,inputTEList(klist).name);
                end
                
                TE = readEchoTime(inputTEListCell);
                
                % dcm2niix: B0 is available with tag 'MagneticFieldStrength'
                B0 = get_value_from_file(inputTEListCell{1},'MagneticFieldStrength');
                if ~isempty(B0)
                    B0 = str2double(B0);
                end
                
            else
                inputTEList = dir([input '/*.mat*']);
                if ~isempty(inputTEList)
                    TE = readEchoTime(inputTEList(1));
                else
                    error('No TE file is identified.');
                end
            end
        end
    else
        %% Option 3: input is a directory containing DICOM files
        sepia_addpath('MEDI');
        
        [~,voxelSize,matrixSize,CF,delta_TE,TE,B0_dir]=Read_DICOM(input);
        
        B0 = CF/(gyro*1e6);
    end
end

%% check user defined input and set default
% B0
if isfield(userDefine,'B0')
    if ~isempty(userDefine.B0)
        B0 = userDefine.B0;
    end
end
if isempty(B0)
    warning('No field strength is found. Assuming B0 = 3.');
    B0 = 3;
end

% B0 direction
if isfield(userDefine,'B0_dir')
    if ~isempty(userDefine.B0_dir)
        B0_dir = userDefine.B0_dir;
    end
end
if isempty(B0_dir)
    warning('No main field direction is found. Assuming B0_dir = [0,0,1].');
    B0_dir = [0,0,1];
end

% voxel size
if isfield(userDefine,'voxelSize')
    if ~isempty(userDefine.voxelSize)
        voxelSize = userDefine.voxelSize;
    end
end
if isempty(voxelSize)
    warning('No spatial resolution is found. Assuming voxelSize = [1,1,1].');
    voxelSize = [1,1,1];
end

% TE
if isfield(userDefine,'TE')
    if ~isempty(userDefine.TE)
        TE = userDefine.TE;
    end
end
if isempty(TE) || isnan(TE(1))
    error('No echo time is found. Please provide the TE information.');
end

% central frequency
if isempty(CF)
    CF = B0 * gyro * 1e6;
end

% get echo spacing
if isempty(delta_TE)
    if length(TE) > 1
        % multi-echo
        delta_TE = TE(2) - TE(1);
    else
        % single echo
        delta_TE = TE(1);
    end
end

%% save header
save([outputDir filesep prefix 'header.mat'],'voxelSize','matrixSize','CF','delta_TE',...
        'TE','B0_dir','B0');

end

%% get echo time from files
function TE = readEchoTime(teFullName)
        
    % get extension from input file
    [~,~,ext] = fileparts(teFullName{1});

    % extract TE based on input data extension
    switch lower(ext)
        case '.mat'
            % if mat file then try to load 'TE' directly
            try load(teFullName{1},'TE');  catch; error('No variable named ''TE''.'); end
        case '.txt'
            % if text file the try to read the TEs line by line
            TE = readTE_MRIConvert_Text(teFullName{1});
        case '.json'
            % JSON file(s)
            TE = readTE_JSON(teFullName);
    end
    TE = TE(:).';
end

%% get value from file
function target_value = get_value_from_file(fname,target_string)
target_value = [];

% read mode
fid = fopen(fname,'r');
% read file line by line
line = fgetl(fid);
    
% start reading lines
while ischar(line)
    % look for target string
    if ~isempty(strfind(line, target_string))
        target_value = get_str(line);
        break
    end
    % start the next line
    line = fgetl(fid);
end

fclose(fid);
end

%% Get value of the tag 
function str=get_str(list_info)
    % find chars ': '
    k_b = strfind(list_info,': ');
    % account for two chars ':' and ' '
    str=(list_info(k_b(1)+2:end));
end
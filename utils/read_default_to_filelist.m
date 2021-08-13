%% inputNIFTIList = read_default_to_filelist(inputDir, filePattern)
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
% Date created: 11 August 2021
% Date modified:
%
%
function [isLoadSuccessful, inputNIFTIList] = read_default_to_filelist(inputDir, filePattern)

disp('#####################################################');
disp('# Checking input directory for SEPIA default format #');
disp('#####################################################');

% check and get filenames
inputNIFTIList = struct();
% filePattern = {'ph','mag','weights','header'}; % don't change the order
isLoadSuccessful = false;

for  k = 1:length(filePattern)
    % get filename
    if ~isempty(filePattern{k})
        if k ~= 4   % NIFTI image input
            [file,numFiles] = get_filename_in_directory(inputDir,filePattern{k},'.nii');
        else        % SEPIA header
            [file,numFiles] = get_filename_in_directory(inputDir,filePattern{k},'.mat');
        end

        % actions given the number of files detected
        if numFiles == 1        % only one file -> get the name

            fprintf('One ''%s'' file is found: %s\n',filePattern{k},file.name);
            inputNIFTIList(k).name = file.name;

        elseif numFiles == 0     % no file -> error

            if k == 1 || k ==4 % essential files: always the 1st and 4th fields 
                warning(['No file with name containing string ''' filePattern{k} ''' is detected.']);
                disp('Reading input directory based on default format failed.');

                return

            else
                disp(['No file with name containing string ''' filePattern{k} ''' is detected.']);
            end

        else % multiple files -> fatal error

            warning(['Multiple files with name containing string ''' filePattern{k} ''' are detected. Make sure the input directory should contain only one file with string ''' filePattern{k} '''.']);
            disp('Reading input directory based on default format failed.');

            return

        end
    end
end

% if no warning reached, then loading is successful
isLoadSuccessful = true;
    
end
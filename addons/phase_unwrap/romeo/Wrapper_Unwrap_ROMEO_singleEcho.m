%% [unwrappedField] = Wrapper_Unwrap_ROMEO_singleEcho(wrappedField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
%
% Input
% --------------
% totalField    : total field map (background + tissue fields), in Hz
% mask          : signal mask
% matrixSize    : size of the input image
% voxelSize     : spatial resolution of each dimension of the data, in mm
% algorParam    : structure contains fields with algorithm-specific parameter(s)
% headerAndExtraData : structure contains extra header info/data for the algorithm
%
% Output
% --------------
% RDF           : local field map
%
% Description: This is a wrapper function to access ROMEO single-echo unwrapping for SEPIA
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 8 March 2020 (v0.8.0)
% Date modified: 16 August 2021 (v1.0, KC)
%
%
function [unwrappedField] = Wrapper_Unwrap_ROMEO_singleEcho(wrappedField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
sepia_universal_variables;
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);

% add path
sepia_addpath('MRITOOLS');

% get algorithm parameters
parameters = check_and_set_algorithm_default(headerAndExtraData, mask);

%% create tmp directory
% Should create a suitable temporary directory on every machine
% Suggestion 20220912 KC: can use the output directory as temporary
% directory for ROMEO output, this should work for SEPIA v1.1
if isfield(headerAndExtraData, 'outputDirectory')
    parameters.output_dir = fullfile(headerAndExtraData.outputDirectory,'romeo_tmp');
else
    parameters.output_dir = fullfile(tempdir, 'romeo_tmp');
end
mkdir(parameters.output_dir);

%% check shared library
% 20231009: KC for latest Matlab versions
if ~ispc
    filepath = fileparts( which('ROMEO.m'));
    libunwind_file = fullfile(filepath,'..','lib','julia','libunwind.so');
    path = getenv('LD_PRELOAD');
    if ~contains(path,libunwind_file)
        setenv('LD_PRELOAD',strcat(path,':',libunwind_file));
    end
end

%% main
[unwrappedField, ~] = ROMEO(wrappedField, parameters);

%% Delete tmp directory
% Remove all temp output files and the temp folder
rmdir(parameters.output_dir, 's')
       
end

%% set parameters
function parameters = check_and_set_algorithm_default(headerAndExtraData, mask)

magn = get_variable_from_headerAndExtraData(headerAndExtraData, 'magnitude', size(mask));

parameters.TE = false;
parameters.no_unwrapped_output = false;
parameters.calculate_B0 = false;
parameters.mag = magn(:,:,:,1); % headerAndExtraData.magn(:,:,:,1); % use first magnitude
parameters.phase_offset_correction = 'off';
parameters.mask = mask;

end

%% ext = get_nifti_extension_from_input(input)
%
% Input
% --------------
% input     : input directory (char/string) or structure containing input
%             NIfTI filenames, as accepted by the SEPIA I/O wrappers
%
% Output
% --------------
% ext       : '.nii' if the input data uses uncompressed NIfTI,
%             otherwise '.nii.gz' (also the default when no NIfTI file
%             can be found, e.g. input is empty)
%
% Description: Detects whether SEPIA outputs should be saved as '.nii' or
%              '.nii.gz' by inspecting the input data, instead of always
%              forcing '.nii.gz'.
%
% Kwok-shing Chan @ DCCN
% kwokshing.chan@donders.ru.nl
% Date created: 23 August 2026 (v1.3.0)
%
function ext = get_nifti_extension_from_input(input)

% default
ext = '.nii.gz';

if isstruct(input)

    for k = 1:numel(input)
        if ~isempty(input(k).name)
            ext = get_extension_from_filename(input(k).name, ext);
            return
        end
    end

elseif (ischar(input) || isstring(input)) && ~isempty(input)

    input = char(input);

    if exist(input,'dir') == 7
        % input is a directory: look for the first NIfTI file inside
        niiList = [dir(fullfile(input,'*.nii.gz')); dir(fullfile(input,'*.nii'))];
        if ~isempty(niiList)
            ext = get_extension_from_filename(niiList(1).name, ext);
        end
    elseif exist(input,'file') == 2
        ext = get_extension_from_filename(input, ext);
    end

end

end

%% get '.nii' or '.nii.gz' from a filename, keeping 'defaultExt' if neither matches
function ext = get_extension_from_filename(filename, defaultExt)

if endsWith(filename, '.nii.gz')
    ext = '.nii.gz';
elseif endsWith(filename, '.nii')
    ext = '.nii';
else
    ext = defaultExt;
end

end

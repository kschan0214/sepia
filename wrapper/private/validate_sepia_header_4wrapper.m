%% Standard script to validate input sepia header
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 22 Jan 2021
% Date modified: 12 August 2021 (v1.0)
%
%
function sepia_header = validate_sepia_header_4wrapper(sepia_header, outputNiftiTemplate)

sepia_universal_variables;

disp('Validating SEPIA header...')

% Magnetic field strength and TE are rarely stored in NIfTI
% Without these two variables -> fatal error
if ~isfield(sepia_header, 'B0')
    error('Variable ''B0'' is missing in the SEPIA header. Please check the input header .mat file.');
end
if ~isfield(sepia_header, 'TE')
    error('Variable ''TE'' is missing in the SEPIA header. Please check the input header .mat file.');
end

% In case some parameters are missing in the header file try to get as
% much information as possible from NIfTI
if ~isfield(sepia_header, 'matrixSize')
    disp('Variable ''matrixSize'' is missing in the SEPIA header. NIfTI header info will be used.')
    sepia_header.matrixSize = outputNiftiTemplate.hdr.dime.dim(2:4);
end
if ~isfield(sepia_header, 'voxelSize')
    disp('Variable ''voxelSize'' is missing in the SEPIA header. NIfTI header info will be used.')
    sepia_header.voxelSize = outputNiftiTemplate.hdr.dime.pixdim(2:4);
end
if ~isfield(sepia_header, 'B0_dir')
    disp('Variable ''B0_dir'' is missing in the SEPIA header. NIfTI header info will be used.')
    sepia_header.B0_dir = get_B0_dir_from_nifti(outputNiftiTemplate);
end
if ~isfield(sepia_header, 'CF')
    sepia_header.CF = sepia_header.B0 * gyro * 1e6;
end
if ~isfield(sepia_header,'delta_TE')
    if length(sepia_header.TE) == 1
        sepia_header.delta_TE = sepia_header.TE;
    else
        sepia_header.delta_TE = sepia_header.TE(2) - sepia_header.TE(1);
    end
end

% make sure B0_dir is a unit vector
sepia_header.B0_dir      = sepia_header.B0_dir ./ norm(sepia_header.B0_dir);
sepia_header.matrixSize  = sepia_header.matrixSize(:).';          % row vectors
sepia_header.voxelSize   = sepia_header.voxelSize(:).';           % row vectors

disp('Input SEPIA header is valid.')

end
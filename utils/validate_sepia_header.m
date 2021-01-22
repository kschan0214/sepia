%% Standard script to validate input sepia header

disp('Validating SEPIA header...')

% Magnetic field strength and TE are rarely stored in NIfTI
% Without these two variables -> fatal error
if ~exist('B0','var')
    error('Variable ''B0'' is missing in the SEPIA header. Please check the input header .mat file.');
end
if ~exist('TE','var')
    error('Variable ''TE'' is missing in the SEPIA header. Please check the input header .mat file.');
end

% In case some parameters are missing in the header file try to get as
% much information as possible from NIfTI
if ~exist('matrixSize','var')
    matrixSize = outputNiftiTemplate.hdr.dime.dim(2:4);
end
if ~exist('voxelSize','var')
    voxelSize = outputNiftiTemplate.hdr.dime.pixdim(2:4);
end
if ~exist('B0_dir','var')
    B0_dir = get_B0_dir_from_nifti(outputNiftiTemplate);
end
if ~exist('CF','var')
    CF = B0 * gyro * 1e6;
end
if ~exist('delta_TE','var')
    if length(TE) ==1
        delta_TE = TE;
    else
        delta_TE = TE(2) - TE(1);
    end
end

% make sure B0_dir is a unit vector
B0_dir      = B0_dir ./ norm(B0_dir);
matrixSize  = matrixSize(:).';          % row vectors
voxelSize   = voxelSize(:).';           % row vectors

disp('Input SEPIA header is valid.')

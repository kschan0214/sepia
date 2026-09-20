%% filterSpec = sepia_get_nifti_uigetfile_filter()
%
% Output
% --------------
% filterSpec    : filter spec to be passed as the 2nd input of uigetfile
%                 for selecting a NIfTI (*.nii or *.nii.gz) file
%
% Description: macOS's native uigetfile file-type selector does not
%              reliably honour a multi-row filter cell array here - it can
%              default to a row that filters out *.nii.gz files, showing
%              them as greyed out/unselectable regardless of row order.
%              A single plain string filter has no such type selector, so
%              it is used on macOS; Windows/Linux keep the more
%              descriptive multi-row filter.
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 19 September 2026
%
%
function filterSpec = sepia_get_nifti_uigetfile_filter()

if ismac
    filterSpec = '*.*';
else
    filterSpec = {'*.nii;*.nii.gz','NIfTI file (*.nii,*.nii.gz)'; ...
                   '*.nii','NIfTI file (*.nii)'; ...
                   '*.nii.gz','Gzipped NIfTI file (*.nii.gz)'};
end

end

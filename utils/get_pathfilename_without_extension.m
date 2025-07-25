%% name = getFilenameWithoutExtension(filepath)
% 
% Input
% --------------
% filepath      : file path (.nii or .nii.gz)
%
% Output
% --------------
% name          : filename without any extension 
%
% Description: get fielname without extension
%
% Kwok-shing Chan @ MGH
% kchan2@mgh.harvard.edu
% Date created: 14 July 2025
% Date modified: 
%
%
function name = get_pathfilename_without_extension(filepath)
[path, name, ext] = fileparts(filepath);
    if strcmp(ext, '.gz')
        % Handle .nii.gz
        [~, name, ~] = fileparts(name); % remove .nii from .nii.gz
    end
    name = fullfile(path,name);
end

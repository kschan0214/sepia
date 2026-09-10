function setup_sepia()
% One-time setup: creates local, machine-specific copies of SEPIA's
% config files from their tracked templates, for whichever ones don't
% already exist.

thisDir = fileparts(mfilename('fullpath'));

% {template filename, local filename, message shown when it's created}
configs = { ...
    'SpecifyToolboxesDirectory.template.m', 'SpecifyToolboxesDirectory.m', ...
        'Please edit this file to point to your local MEDI/STI Suite/FANSI/SEGUE installations,\nthen re-run sepia_addpath.'; ...
    'SpecifyAtlasDirectory.template.m',     'SpecifyAtlasDirectory.m', ...
        'Please edit this file if your atlas folders are not at download_atlas.m''s default locations.' ...
    };

for k = 1:size(configs,1)
    templateFile = fullfile(thisDir, configs{k,1});
    localFile    = fullfile(thisDir, configs{k,2});

    if ~isfile(localFile)
        copyfile(templateFile, localFile);
        fprintf(['SEPIA setup: created %s from template.\n' configs{k,3} '\n'], configs{k,2});
%         edit(localFile);   % opens it in the MATLAB editor for convenience
%     else
%         fprintf('%s already exists, leaving it untouched.\n', configs{k,2});
    end
end

end

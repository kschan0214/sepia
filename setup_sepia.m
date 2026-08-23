function setup_sepia()
% One-time setup: creates a local, machine-specific copy of the
% toolbox path configuration file if it doesn't already exist.

thisDir     = fileparts(mfilename('fullpath'));
templateFile = fullfile(thisDir, 'SpecifyToolboxesDirectory.template.m');
localFile    = fullfile(thisDir, 'SpecifyToolboxesDirectory.m');

if ~isfile(localFile)
    copyfile(templateFile, localFile);
    fprintf(['SEPIA setup: created SpecifyToolboxesDirectory.m from template.\n' ...
              'Please edit this file to point to your local MEDI/STI Suite/FANSI/SEGUE installations,\n' ...
              'then re-run sepia_addpath.\n']);
%     edit(localFile);   % opens it in the MATLAB editor for convenience
% else
%     fprintf('SpecifyToolboxesDirectory.m already exists, leaving it untouched.\n');
end

end
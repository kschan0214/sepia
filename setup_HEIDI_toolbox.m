function setup_HEIDI_toolbox()
% One-time setup: check whether the HEIDI_SEPIAready package is available;
% if not, download it into SEPIA_HOME/external/ and register its path in
% SpecifyToolboxesDirectory.m.
%
% This script only runs when you call it explicitly - it is not invoked
% automatically by sepia_addpath, nor by the HEIDI/LSQR+HEIDI wrappers
% themselves (which just read whatever HEIDI_HOME is already configured).
%
% Kwok-shing Chan @ MGH
% kchan2@mgh.harvard.edu
% Date created: 6 September 2026
%
%
SEPIA_HOME = fileparts(mfilename('fullpath'));
addpath(SEPIA_HOME);

% TODO: fill in once HEIDI_SEPIAready is uploaded somewhere (e.g. a GitHub
% Release asset on the SEPIA repo, or a Zenodo record) - see the archive
% expected to unzip directly into a folder containing 'HEIDI/' and
% 'LSQR/' subfolders (i.e. matching the existing external/HEIDI_SEPIAready
% sibling-folder convention).
url = '<HEIDI_SEPIAready_DOWNLOAD_URL>';

% make sure SpecifyToolboxesDirectory.m exists
setup_sepia();

% get the currently configured HEIDI_HOME, if any
run(fullfile(SEPIA_HOME,'SpecifyToolboxesDirectory.m'));
if ~exist('HEIDI_HOME','var')
    HEIDI_HOME = [];
end

% marker file used to confirm a valid HEIDI_SEPIAready installation
markerFile = fullfile('HEIDI','GradientAnisotropicDiffusionImageFilter');

if ~isempty(HEIDI_HOME) && exist(fullfile(HEIDI_HOME,markerFile),'file') == 2
    fprintf('HEIDI_SEPIAready is already available at %s\n', HEIDI_HOME);
    return
end

if isempty(url) || strcmp(url,'<HEIDI_SEPIAready_DOWNLOAD_URL>')
    error('SEPIA:HEIDIDownloadURLNotSet', ...
        'The download URL for HEIDI_SEPIAready has not been set in setup_HEIDI_toolbox.m yet. Please obtain the package manually and set HEIDI_HOME in SpecifyToolboxesDirectory.m (or the Utility tab''s Manage Dependency panel) instead.');
end

fprintf('HEIDI_SEPIAready not found. Downloading from %s...\n', url);

downloadDir = fullfile(SEPIA_HOME,'external','HEIDI_SEPIAready');
if exist(downloadDir,'dir') ~= 7
    mkdir(downloadDir);
end

zipFile = fullfile(SEPIA_HOME,'external','HEIDI_SEPIAready.zip');

status = system(sprintf('wget -O %s --no-check-certificate %s', zipFile, url));
if status ~= 0 || exist(zipFile,'file') ~= 2
    error('SEPIA:HEIDIDownloadFailed', ...
        'Failed to download HEIDI_SEPIAready from %s. Please check your network connection or download it manually.', url);
end

status = system(sprintf('unzip -o %s -d %s', zipFile, downloadDir));
delete(zipFile);

HEIDI_HOME = downloadDir;

if status ~= 0 || exist(fullfile(HEIDI_HOME,markerFile),'file') ~= 2
    error('SEPIA:HEIDIExtractFailed', ...
        'Failed to extract HEIDI_SEPIAready to %s. Please check the download or extract it manually.', downloadDir);
end

% grant execute permission to the bundled binary
system(sprintf('chmod a+x %s', fullfile(HEIDI_HOME,markerFile)));

fprintf('HEIDI_SEPIAready downloaded to %s\n', HEIDI_HOME);

% register the new path in SpecifyToolboxesDirectory.m
isUpdated = update_toolbox_directory_entry(SEPIA_HOME, 'HEIDI_HOME', HEIDI_HOME);
if isUpdated
    fprintf('SpecifyToolboxesDirectory.m updated with the new HEIDI_HOME.\n');
else
    fprintf('SpecifyToolboxesDirectory.m already points to this HEIDI_HOME.\n');
end

end

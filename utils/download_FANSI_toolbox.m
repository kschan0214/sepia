function download_FANSI_toolbox()
% One-time setup: check whether FANSI-toolbox (https://gitlab.com/cmilovic/FANSI-toolbox)
% is available; if not, download the pinned commit into SEPIA_HOME/external/
% and register its path in SpecifyToolboxesDirectory.m.
%
% This script only runs when you call it explicitly - it is not invoked
% automatically by sepia_addpath.
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 23 August 2026
%
%
SEPIA_HOME = fileparts(fileparts(mfilename('fullpath')));
addpath(SEPIA_HOME);

% pinned commit of FANSI-toolbox
COMMIT       = '93044a301e2cbe94177c76e8add77c9859c61283';
COMMIT_SHORT = COMMIT(1:8);

% make sure SpecifyToolboxesDirectory.m exists
setup_sepia();

% get the currently configured FANSI_HOME, if any
run(fullfile(SEPIA_HOME,'SpecifyToolboxesDirectory.m'));
if ~exist('FANSI_HOME','var')
    FANSI_HOME = [];
end

% marker file used to confirm a valid FANSI-toolbox installation
markerFile = 'dipole_kernel_fansi.m';

if ~isempty(FANSI_HOME) && exist(fullfile(FANSI_HOME,markerFile),'file') == 2
    fprintf('FANSI-toolbox is already available at %s\n', FANSI_HOME);
    return
end

fprintf('FANSI-toolbox not found. Downloading commit %s from GitLab...\n', COMMIT_SHORT);

downloadDir = fullfile(SEPIA_HOME,'external','FANSI_toolbox');
if exist(downloadDir,'dir') ~= 7
    mkdir(downloadDir);
end

zipFile = fullfile(downloadDir, ['FANSI-toolbox-' COMMIT_SHORT '.zip']);
url     = sprintf('https://gitlab.com/cmilovic/FANSI-toolbox/-/archive/%s/FANSI-toolbox-%s.zip', COMMIT, COMMIT_SHORT);

status = system(sprintf('wget -O %s --no-check-certificate %s', zipFile, url));
if status ~= 0 || exist(zipFile,'file') ~= 2
    error('SEPIA:FANSIDownloadFailed', ...
        'Failed to download FANSI-toolbox from %s. Please check your network connection or download it manually.', url);
end

status = system(sprintf('unzip -o %s -d %s', zipFile, downloadDir));
delete(zipFile);

FANSI_HOME = fullfile(downloadDir, ['FANSI-toolbox-' COMMIT_SHORT], filesep);

if status ~= 0 || exist(fullfile(FANSI_HOME,markerFile),'file') ~= 2
    error('SEPIA:FANSIExtractFailed', ...
        'Failed to extract FANSI-toolbox to %s. Please check the download or extract it manually.', downloadDir);
end

fprintf('FANSI-toolbox downloaded to %s\n', FANSI_HOME);

% register the new path in SpecifyToolboxesDirectory.m
isUpdated = update_toolbox_directory_entry(SEPIA_HOME, 'FANSI_HOME', FANSI_HOME);
if isUpdated
    fprintf('SpecifyToolboxesDirectory.m updated with the new FANSI_HOME.\n');
else
    fprintf('SpecifyToolboxesDirectory.m already points to this FANSI_HOME.\n');
end

end

function setup_tMPPCA_toolbox()
% Check whether Tensor-MP-PCA (https://github.com/Neurophysics-CFIN/Tensor-MP-PCA)
% is available on the path; if not, download it into SEPIA_HOME/external/
% and add it to the path.
%
% Unlike FANSI/HEIDI, Tensor-MP-PCA needs no SpecifyToolboxesDirectory.m
% entry - its location is always SEPIA_HOME/external/Tensor-MP-PCA, so
% nothing is registered there; this function just makes sure the code is
% present and on the path.
%
% This function is called both automatically (lazily, from
% SepiaIOWrapper.m/UnwrapPhaseMacroIOWrapper.m, the first time denoising is
% actually requested) and can be called directly by the user, e.g. via
% setup_sepia_downloads() to pre-fetch every auto-downloadable toolbox at
% once.
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 September 2026
%
%
SEPIA_HOME = fileparts(fileparts(mfilename('fullpath')));

tMPPCA_HOME = fullfile(SEPIA_HOME,'external','Tensor-MP-PCA');
if exist(tMPPCA_HOME,'dir'); addpath(genpath(tMPPCA_HOME)); end

if exist('denoise_recursive_tensor', 'file')
    fprintf('Tensor-MP-PCA is already available at %s\n', tMPPCA_HOME);
    return
end

fprintf('Tensor-MP-PCA not found. Downloading from GitHub...\n');

if exist(tMPPCA_HOME,'dir') ~= 7
    mkdir(tMPPCA_HOME);
end

zipFile = strcat(tMPPCA_HOME,'.zip');
url     = 'https://github.com/Neurophysics-CFIN/Tensor-MP-PCA/archive/refs/heads/main.zip';

status = system(sprintf('wget -O %s --no-check-certificate %s', zipFile, url));
if status ~= 0 || exist(zipFile,'file') ~= 2
    error('SEPIA:tMPPCADownloadFailed', ...
        'Failed to download Tensor-MP-PCA from %s. Please check your network connection or download it manually.', url);
end

status = system(sprintf('unzip -o %s -d %s', zipFile, tMPPCA_HOME));
delete(zipFile);

addpath(genpath(tMPPCA_HOME));

if status ~= 0 || ~exist('denoise_recursive_tensor', 'file')
    error('SEPIA:tMPPCAExtractFailed', ...
        'Failed to extract Tensor-MP-PCA to %s. Please check the download or extract it manually.', tMPPCA_HOME);
end

fprintf('Tensor-MP-PCA downloaded to %s\n', tMPPCA_HOME);

end

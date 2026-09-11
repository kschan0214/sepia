function download_atlas()
% Cross-platform (Windows/macOS/Linux) atlas downloader - replaces
% download_atlas.sh, which cannot run on Windows at all (no bash) and is
% unreliable on macOS (its `readlink -f` call depends on a GNU coreutils
% flag BSD/macOS readlink doesn't support).
%
% Downloads the CIT168, MuSus-100 and AHEAD atlases into SEPIA_HOME/atlas/,
% in the exact layout SpecifyAtlasDirectory.m expects
% (CIT168_Reinf_Learn_v1.1.0/, MuSus-100_Atlas/, AHEAD_atlas/). Uses only
% MATLAB's built-in websave/unzip/untar - no curl/wget/tar/git command-line
% tools required, so it works identically on every platform MATLAB itself
% runs on.
%
% Safe to run repeatedly: each atlas is skipped if its destination folder
% already exists and isn't empty. A failure on one atlas (e.g. no network,
% or a source going offline) does not stop the others from being
% attempted; a summary of what succeeded/failed is printed at the end.
%
% Kwok-shing Chan
% Date created: 10 September 2026
%
SEPIA_HOME = fileparts(mfilename('fullpath'));
atlasDir = fullfile(SEPIA_HOME, 'atlas');
if ~isfolder(atlasDir); mkdir(atlasDir); end

atlases  = {'CIT168', 'MuSus-100', 'AHEAD'};
downloadFcns = {@download_CIT168, @download_MuSus100, @download_AHEAD};

isOK = false(size(atlases));
for k = 1:numel(atlases)
    fprintf('\n=== %s atlas ===\n', atlases{k});
    try
        downloadFcns{k}(atlasDir);
        isOK(k) = true;
    catch ME
        fprintf(2, 'Failed to download %s atlas: %s\n', atlases{k}, ME.message);
    end
end

fprintf('\n=== Summary ===\n');
for k = 1:numel(atlases)
    if isOK(k)
        fprintf('%-12s OK\n', atlases{k});
    else
        fprintf('%-12s FAILED (see message above)\n', atlases{k});
    end
end

end

%% CIT168 reinforcement-learning atlas (OSF)
function download_CIT168(atlasDir)

destDir = fullfile(atlasDir, 'CIT168_Reinf_Learn_v1.1.0');
if isfolder(destDir) && ~isempty(dir(destDir))
    fprintf('CIT168 atlas already present at %s\n', destDir);
    return
end
mkdir(destDir);

url = 'https://files.osf.io/v1/resources/jkzwp/providers/osfstorage/5b11f8d6f1f288000d6343aa/?zip=';
zipFile = fullfile(atlasDir, 'CIT168_Reinf_Learn_v1.1.0.zip');

fprintf('Downloading CIT168 atlas...\n');
websave(zipFile, url, weboptions('CertificateFilename', ''));
unzip(zipFile, destDir);
delete(zipFile);

fprintf('CIT168 atlas downloaded to %s\n', destDir);

end

%% MuSus-100 atlas (GitHub - downloaded as a zip archive, not `git clone`,
%% so this needs no git installation)
function download_MuSus100(atlasDir)

destDir = fullfile(atlasDir, 'MuSus-100_Atlas');
if isfolder(destDir) && ~isempty(dir(destDir))
    fprintf('MuSus-100 atlas already present at %s\n', destDir);
    return
end

url = 'https://github.com/SMILE-Lab-ShanghaiTech/MuSus-100_Atlas/archive/refs/heads/main.zip';
zipFile = fullfile(atlasDir, 'MuSus-100_Atlas.zip');

fprintf('Downloading MuSus-100 atlas...\n');
websave(zipFile, url, weboptions('CertificateFilename', ''));

tmpDir = fullfile(atlasDir, 'MuSus-100_Atlas_tmp');
unzip(zipFile, tmpDir);
delete(zipFile);

% a GitHub archive zip extracts into a single '<repo>-<branch>' subfolder;
% move its contents up to the expected destDir
extracted = dir(fullfile(tmpDir, 'MuSus-100_Atlas-*'));
extracted = extracted([extracted.isdir]);
if isempty(extracted)
    rmdir(tmpDir, 's');
    error('SEPIA:MuSus100UnexpectedLayout', ...
        'Unexpected zip layout after extracting the MuSus-100 atlas archive.');
end
movefile(fullfile(tmpDir, extracted(1).name), destDir);
rmdir(tmpDir, 's');

fprintf('MuSus-100 atlas downloaded to %s\n', destDir);

end

%% AHEAD atlas (figshare - two tar.gz archives extracted into the same folder)
function download_AHEAD(atlasDir)

destDir = fullfile(atlasDir, 'AHEAD_atlas');
if isfolder(destDir) && ~isempty(dir(destDir))
    fprintf('AHEAD atlas already present at %s\n', destDir);
    return
end
mkdir(destDir);

files = { ...
    'structure_mni09b.tar.gz', 'https://uvaauas.figshare.com/ndownloader/files/21209229'; ...
    'Templates_mni09b.tar.gz', 'https://uvaauas.figshare.com/ndownloader/files/21209235' ...
    };

for i = 1:size(files,1)
    tarFile = fullfile(atlasDir, files{i,1});
    fprintf('Downloading %s...\n', files{i,1});
    websave(tarFile, files{i,2}, weboptions('CertificateFilename', ''));
    untar(tarFile, destDir);
    delete(tarFile);
end

fprintf('AHEAD atlas downloaded to %s\n', destDir);

end

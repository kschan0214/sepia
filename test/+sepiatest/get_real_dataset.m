%% ds = sepiatest.get_real_dataset()
%
% Description: locates the real dataset used by the Tier-2 toolbox-dependent
% regression matrix. Looked up from (in order):
%   1. environment variables SEPIA_TEST_REAL_DATA_DIR (input directory) and
%      SEPIA_TEST_REAL_DATA_MASK (brain mask NIfTI)
%   2. a test/config/real_dataset.json file (not tracked in git - see
%      test/.gitignore), with fields "inputDir" and "maskFile", as a
%      convenience so a developer doesn't have to export env vars every
%      shell session
%
% Returns ds.inputDir='' / ds.maskFile='' (both empty) if not configured,
% or if the configured paths don't actually exist - callers should treat
% that as "Tier 2 not available on this machine" and skip gracefully
% (assumeTrue), never fail.
%
% inputDir is passed directly as sepiaIO's `input` (a directory containing
% phase/magnitude NIfTI+JSON files, BIDS-style or SEPIA's own naming
% convention - see docs/getting_started/Data-preparation.rst); maskFile is
% passed as sepiaIO's `maskFullName` and need not live inside inputDir.
%
function ds = get_real_dataset()

ds = struct('inputDir', '', 'maskFile', '');

ds.inputDir  = getenv('SEPIA_TEST_REAL_DATA_DIR');
ds.maskFile  = getenv('SEPIA_TEST_REAL_DATA_MASK');

if isempty(ds.inputDir) || isempty(ds.maskFile)
    cfgFile = fullfile(sepiatest.test_root(), 'config', 'real_dataset.json');
    if isfile(cfgFile)
        cfg = jsondecode(fileread(cfgFile));
        if isempty(ds.inputDir) && isfield(cfg,'inputDir')
            ds.inputDir = cfg.inputDir;
        end
        if isempty(ds.maskFile) && isfield(cfg,'maskFile')
            ds.maskFile = cfg.maskFile;
        end
    end
end

if ~isempty(ds.inputDir) && exist(ds.inputDir,'dir') ~= 7
    warning('sepiatest:realDataNotFound', 'Configured real-dataset inputDir does not exist: %s', ds.inputDir);
    ds.inputDir = '';
end
if ~isempty(ds.maskFile) && exist(ds.maskFile,'file') ~= 2
    warning('sepiatest:realDataNotFound', 'Configured real-dataset maskFile does not exist: %s', ds.maskFile);
    ds.maskFile = '';
end

end

function regenerate_reference(varargin)
%% regenerate_reference('confirm', true, 'tier', 'tier1', 'filter', '*')
%
% Description: (re)generates saved reference outputs used by the
% regression comparisons in test/tier1_smoke and test/tier2_matrix.
% This is a DELIBERATE, DESTRUCTIVE action - it overwrites files under
% test/references/ - so it refuses to run at all unless called with
% 'confirm', true, and prints an old-vs-new diff for anything it
% overwrites so the change is reviewable (e.g. via `git diff`) before
% being committed.
%
% Name-value arguments
% --------------
% confirm : (required) must be passed as true, else errors immediately.
% tier    : 'tier1' (default) or 'tier2'. Only 'tier1' (the toolbox-free
%           smoke tier) is implemented so far.
% filter  : glob-style filter on which reference file(s) to regenerate
%           within the tier, e.g. 'smoke_TKD' (default '*' = all).
%

p = inputParser;
p.addParameter('confirm', false);
p.addParameter('tier', 'tier1');
p.addParameter('filter', '*');
p.parse(varargin{:});
opt = p.Results;

if ~isequal(opt.confirm, true)
    error('regenerate_reference:confirmRequired', ...
        ['Refusing to overwrite reference outputs without ''confirm'',true. ', ...
         'This is a deliberate, destructive action - re-run as ', ...
         'regenerate_reference(''confirm'',true) only after reviewing why the ', ...
         'existing reference is expected to change.']);
end

testRoot = sepiatest.test_root();
SEPIA_HOME = sepiatest.sepia_home();
addpath(SEPIA_HOME);
sepia_addpath;
addpath(testRoot);
addpath(fullfile(testRoot,'phantom'));

switch lower(opt.tier)
    case 'tier1'
        regenerate_tier1(testRoot, opt.filter);
    case 'tier2'
        regenerate_tier2(testRoot, opt.filter);
    otherwise
        error('regenerate_reference:unsupportedTier', ...
            'Tier "%s" is not implemented yet in regenerate_reference.m.', opt.tier);
end

end

%% Tier 1: toolbox-free smoke reference (mirrors TestSmokePhantom.m exactly)
function regenerate_tier1(testRoot, filterExpr)

qsmMethods = {'TKD', 'Closed-form solution', 'iLSQR'};

tmp = tempname; mkdir(tmp);
phantomPaths = generate_synthetic_phantom(fullfile(tmp,'phantom'));

refDir = fullfile(testRoot, 'references', 'tier1');
if exist(refDir,'dir') ~= 7; mkdir(refDir); end

for k = 1:numel(qsmMethods)
    qsmMethod = qsmMethods{k};
    slug = matlab.lang.makeValidName(qsmMethod);
    if ~isempty(filterExpr) && ~strcmp(filterExpr,'*') && ~strcmpi(filterExpr, slug)
        continue
    end

    algorParam = struct();
    algorParam.general.isInvert = false;
    algorParam.general.isBET    = false;
    algorParam.unwrap.unwrapMethod   = 'None';
    algorParam.unwrap.echoCombMethod = 'Optimum weights';
    algorParam.bfr.method = 'VSHARP';
    algorParam.qsm.method = qsmMethod;

    outputPrefix = fullfile(tmp, ['out_' slug], 'sepia');
    sepiaIO(phantomPaths.input, outputPrefix, phantomPaths.mask, algorParam);

    % sepiaIO strips test/ off the path internally (see TestSmokePhantom.m
    % comment) - re-add before using sepiatest.* below.
    addpath(testRoot); addpath(fullfile(testRoot,'tools'));

    maskImg = load_nii_img_only(phantomPaths.mask) > 0;
    chiImg  = load_nii_img_only([outputPrefix '_Chimap.nii.gz']);
    stats   = sepiatest.mask_stats(chiImg, maskImg);

    refFile = fullfile(refDir, sprintf('smoke_%s.mat', slug));

    oldStats = [];
    if isfile(refFile)
        old = load(refFile);
        oldStats = old.stats;
    end
    print_diff(qsmMethod, oldStats, stats);

    meta = struct();
    meta.sepiaVersion  = get_sepia_version();
    meta.generatedDate = datestr(datetime('now'),'yyyy-mm-ddTHH:MM:SS');
    meta.matlabVersion = version();
    meta.hostname      = getenv_or('HOSTNAME','');

    save(refFile, 'stats', 'meta');
    fprintf('Wrote %s\n', refFile);
end

end

%% Tier 2: real-dataset method matrices (mirrors TestQSMMatrix.m /
%% TestBFRMatrix.m / TestUnwrapMatrix.m exactly - one category per call)
function regenerate_tier2(testRoot, filterExpr)

ds = sepiatest.get_real_dataset();
if isempty(ds.inputDir) || isempty(ds.maskFile)
    error('regenerate_reference:noRealDataset', ...
        ['No real dataset configured (SEPIA_TEST_REAL_DATA_DIR/SEPIA_TEST_REAL_DATA_MASK or ', ...
         'test/config/real_dataset.json) - cannot regenerate Tier-2 references.']);
end

toolboxes = sepiatest.discover_toolboxes();
addpath(testRoot); addpath(fullfile(testRoot,'tools')); % discover_toolboxes calls sepia_addpath internally, which strips test/ off the path

sepia_universal_variables;

refDir = fullfile(testRoot, 'references', 'tier2');
if exist(refDir,'dir') ~= 7; mkdir(refDir); end

% category, method list, prefix, per-method toolbox-key function, and the
% algorParam field this category varies (holding the other two stages at
% their toolbox-free defaults, same as the corresponding Test*.m class)
categories = { ...
    struct('name','qsm',    'methods', {methodQSMName(:)'},    'prefix','qsm_',    'toolboxKeyFn', @sepiatest.qsm_method_toolbox_key,    'field','qsm.method',    'outputKey','Chimap',    'outputDesc','QSM'), ...
    struct('name','bfr',    'methods', {methodBFRName(:)'},    'prefix','bfr_',    'toolboxKeyFn', @sepiatest.bfr_method_toolbox_key,    'field','bfr.method',    'outputKey','localfield','outputDesc','BFR'), ...
    struct('name','unwrap', 'methods', {methodUnwrapName(:)'}, 'prefix','unwrap_', 'toolboxKeyFn', @sepiatest.unwrap_method_toolbox_key, 'field','unwrap.unwrapMethod', 'outputKey','fieldmap', 'outputDesc','Unwrap') ...
    };

for c = 1:numel(categories)
    cat = categories{c};
    for k = 1:numel(cat.methods)
        method = cat.methods{k};
        slug = matlab.lang.makeValidName(method);
        if ~isempty(filterExpr) && ~strcmp(filterExpr,'*') && ~strcmpi(filterExpr, slug)
            continue
        end

        toolboxKey = cat.toolboxKeyFn(method);
        if strcmp(toolboxKey, 'unavailable')
            fprintf('--- %s:%s --- skipped (no model/checkpoint files available)\n', cat.outputDesc, method);
            continue
        elseif strcmp(toolboxKey, 'excluded')
            fprintf('--- %s:%s --- skipped (deliberately excluded, see test/README.md)\n', cat.outputDesc, method);
            continue
        elseif ~strcmp(toolboxKey, 'none') && ~toolboxes.(toolboxKey)
            fprintf('--- %s:%s --- skipped (%s toolbox not installed on this machine)\n', cat.outputDesc, method, toolboxKey);
            continue
        end

        tmp = tempname; mkdir(tmp);
        outputPrefix = fullfile(tmp, 'sepia');

        algorParam = struct();
        algorParam.general.isInvert = false;
        algorParam.general.isBET    = false;
        algorParam.unwrap.unwrapMethod   = 'None';
        algorParam.unwrap.echoCombMethod = 'Optimum weights';
        algorParam.bfr.method = 'VSHARP';
        algorParam.qsm.method = 'TKD';
        algorParam = setfield_dotted(algorParam, cat.field, method); %#ok<*AGROW>
        if strcmp(cat.name,'qsm') && strcmp(method, 'NDI')
            algorParam.qsm.isGPU = false;
        end

        sepiaIO(ds.inputDir, outputPrefix, ds.maskFile, algorParam);
        addpath(testRoot); addpath(fullfile(testRoot,'tools')); % re-add after sepiaIO's internal sepia_addpath strips it

        maskImg  = load_nii_img_only(ds.maskFile) > 0;
        outImg   = load_nii_img_only([outputPrefix '_' cat.outputKey '.nii.gz']);
        stats    = sepiatest.mask_stats(outImg, maskImg);

        refFile = fullfile(refDir, sprintf('%s%s.mat', cat.prefix, slug));

        oldStats = [];
        if isfile(refFile)
            old = load(refFile);
            oldStats = old.stats;
        end
        print_diff(sprintf('%s:%s', cat.outputDesc, method), oldStats, stats);

        meta = struct();
        meta.sepiaVersion  = get_sepia_version();
        meta.generatedDate = datestr(datetime('now'),'yyyy-mm-ddTHH:MM:SS');
        meta.matlabVersion = version();
        meta.hostname      = getenv_or('HOSTNAME','');
        meta.datasetInputDir = ds.inputDir;
        meta.datasetMaskFile = ds.maskFile;

        save(refFile, 'stats', 'meta');
        fprintf('Wrote %s\n', refFile);
    end
end

end

%% algorParam = setfield_dotted(algorParam, 'qsm.method', 'TKD') sets algorParam.qsm.method='TKD'
function s = setfield_dotted(s, dottedField, value)
parts = strsplit(dottedField, '.');
s.(parts{1}).(parts{2}) = value;
end

function print_diff(label, oldStats, newStats)
fprintf('--- %s ---\n', label);
if isempty(oldStats)
    fprintf('  (no existing reference; creating new)\n');
    return
end
fields = {'mean','std','median','p1','p5','p25','p75','p95','p99'};
for i = 1:numel(fields)
    f = fields{i};
    if isfield(oldStats,f)
        pctChange = 100*(newStats.(f) - oldStats.(f)) / max(abs(oldStats.(f)), eps);
        fprintf('  %-8s old=%.6g  new=%.6g  (%+.3f%%)\n', f, oldStats.(f), newStats.(f), pctChange);
    end
end
end

function v = get_sepia_version()
try
    sepia_universal_variables;
    v = SEPIA_version;
catch
    v = 'unknown';
end
end

function v = getenv_or(name, default)
v = getenv(name);
if isempty(v); v = default; end
end

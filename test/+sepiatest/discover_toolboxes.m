%% tb = sepiatest.discover_toolboxes()
%
% Description: reports which of SEPIA's optional external toolboxes/add-ons
% are actually installed and reachable on this machine, by reading the
% same SpecifyToolboxesDirectory.m that sepia_addpath.m uses (auto-created
% from the template on first call if missing, all paths empty). Used by
% Tier-2 regression tests to skip (not fail) rows whose required toolbox
% isn't present.
%
% Output
% --------------
% tb : struct with logical fields MEDI, STISuite, FANSI, SEGUE, MRITOOLS,
%      MRISC, ANTs, HEIDI, ChiSepNet. HEIDI is also gated on being run on
%      Linux (its wrapper only supports Linux, see docs/method/qsm/HEIDI.rst).
%
function tb = discover_toolboxes()

SEPIA_HOME = sepiatest.sepia_home();

% ensure SpecifyToolboxesDirectory.m exists (auto-created from template by
% setup_sepia() if missing) and load the *_HOME variables it defines,
% without modifying the MATLAB path (isStartCheck=false, method='None').
sepia_addpath('None', false);

% defensively pre-initialise all expected variables to [] in case a
% customised SpecifyToolboxesDirectory.m on some machine doesn't define
% all of them (run() below then just overwrites whichever it does define;
% this function has no nested functions, so introducing/reassigning
% variables here via run() is fine).
MEDI_HOME = []; STISuite_HOME = []; FANSI_HOME = []; SEGUE_HOME = [];
MRITOOLS_HOME = []; MRISC_HOME = []; ANTS_HOME = [];
HEIDI_HOME = []; ChiSepNet_HOME = [];

run(fullfile(SEPIA_HOME,'SpecifyToolboxesDirectory.m'));

% HEIDI_HOME defaults to the historical sibling-folder convention if not
% explicitly configured (kept in sync with the HEIDI wrappers' own fallback).
if isempty(HEIDI_HOME)
    HEIDI_HOME = fullfile(SEPIA_HOME,'..','external','HEIDI_SEPIAready');
end

tb = struct();
tb.MEDI      = ~isempty(MEDI_HOME)     && exist(MEDI_HOME,'dir')     == 7;
tb.STISuite  = ~isempty(STISuite_HOME) && exist(STISuite_HOME,'dir') == 7;
tb.FANSI     = ~isempty(FANSI_HOME)    && exist(FANSI_HOME,'dir')    == 7;
tb.SEGUE     = ~isempty(SEGUE_HOME)    && exist(SEGUE_HOME,'dir')    == 7;
tb.MRITOOLS  = ~isempty(MRITOOLS_HOME) && exist(MRITOOLS_HOME,'dir') == 7;
tb.MRISC     = ~isempty(MRISC_HOME)    && exist(MRISC_HOME,'dir')    == 7;
tb.ANTs      = ~isempty(ANTS_HOME)     && exist(ANTS_HOME,'dir')     == 7;
tb.HEIDI     = isunix && ~ismac && exist(HEIDI_HOME,'dir') == 7;
tb.ChiSepNet = ~isempty(ChiSepNet_HOME) && exist(ChiSepNet_HOME,'dir') == 7;

end

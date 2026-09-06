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
%      MRISC, ANTs, HEIDI (the HEIDI add-on package, checked separately
%      since it is not one of the *_HOME toolboxes but an external folder
%      alongside SEPIA_HOME, see docs/method/qsm/HEIDI.rst)
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

run(fullfile(SEPIA_HOME,'SpecifyToolboxesDirectory.m'));

tb = struct();
tb.MEDI      = ~isempty(MEDI_HOME)     && exist(MEDI_HOME,'dir')     == 7;
tb.STISuite  = ~isempty(STISuite_HOME) && exist(STISuite_HOME,'dir') == 7;
tb.FANSI     = ~isempty(FANSI_HOME)    && exist(FANSI_HOME,'dir')    == 7;
tb.SEGUE     = ~isempty(SEGUE_HOME)    && exist(SEGUE_HOME,'dir')    == 7;
tb.MRITOOLS  = ~isempty(MRITOOLS_HOME) && exist(MRITOOLS_HOME,'dir') == 7;
tb.MRISC     = ~isempty(MRISC_HOME)    && exist(MRISC_HOME,'dir')    == 7;
tb.ANTs      = ~isempty(ANTS_HOME)     && exist(ANTS_HOME,'dir')     == 7;

% HEIDI add-on package lives at <SEPIA_HOME>/../external/HEIDI_SEPIAready
% (hard-coded location, see docs/method/qsm/HEIDI.rst), Linux-only.
tb.HEIDI = isunix && ~ismac && exist(fullfile(SEPIA_HOME,'..','external','HEIDI_SEPIAready'),'dir') == 7;

end

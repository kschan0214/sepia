%% paths = get_sepia_toolbox_home_paths()
%
% Output
% --------------
% paths     : structure with the *_HOME directories declared in
%             SpecifyToolboxesDirectory.m (empty if a given toolbox is
%             not configured there)
%
% Description: Small helper to read the external toolbox directories
%              configured in SpecifyToolboxesDirectory.m, e.g. to check
%              whether a given toolbox is actually available before
%              defaulting the GUI to a method that depends on it.
%
% Kwok-shing Chan @ DCCN
% kwokshing.chan@donders.ru.nl
% Date created: 23 August 2026 (v1.3.0)
%
function paths = get_sepia_toolbox_home_paths()

FANSI_HOME      = [];
MEDI_HOME       = [];
STISuite_HOME   = [];
SEGUE_HOME      = [];
MRITOOLS_HOME   = [];
MRISC_HOME      = [];
ANTS_HOME       = [];

SpecifyToolboxesDirectory;

paths.FANSI_HOME     = FANSI_HOME;
paths.MEDI_HOME      = MEDI_HOME;
paths.STISuite_HOME  = STISuite_HOME;
paths.SEGUE_HOME     = SEGUE_HOME;
paths.MRITOOLS_HOME  = MRITOOLS_HOME;
paths.MRISC_HOME     = MRISC_HOME;
paths.ANTS_HOME      = ANTS_HOME;

end

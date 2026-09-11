function download_toolboxes()
% One-stop setup: checks every SEPIA add-on that can be fetched
% automatically (FANSI, HEIDI, Tensor-MP-PCA) and downloads whichever
% ones are missing, instead of you having to call each
% download_<toolbox>_toolbox() script separately.
%
% Safe to run repeatedly - each individual setup function skips its own
% toolbox if it's already installed. A failure on one toolbox (e.g. no
% network, or a download URL not yet configured) does not stop the others
% from being checked; a summary of what succeeded/failed is printed at
% the end.
%
% Note: this only covers add-ons with a public, scriptable download.
% LPCNN, QSMnet, xQSM, BFRnet and Chi-separation still require manual
% setup (see their own docs pages), since they depend on either a
% model/checkpoint file not publicly hosted, or a toolbox this script has
% no redistribution rights for.
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 September 2026
%
%
SEPIA_HOME = fileparts(mfilename('fullpath'));
addpath(SEPIA_HOME);
addpath(fullfile(SEPIA_HOME,'utils')); % download_FANSI_toolbox/download_HEIDI_toolbox/download_tMPPCA_toolbox live here

toolboxes = {'FANSI', 'HEIDI', 'Tensor-MP-PCA'};
setupFcns = {@download_FANSI_toolbox, @download_HEIDI_toolbox, @download_tMPPCA_toolbox};

isOK = false(size(toolboxes));
for k = 1:numel(toolboxes)
    fprintf('\n=== %s ===\n', toolboxes{k});
    try
        setupFcns{k}();
        isOK(k) = true;
    catch ME
        fprintf(2, 'Failed to set up %s: %s\n', toolboxes{k}, ME.message);
    end
end

fprintf('\n=== Summary ===\n');
for k = 1:numel(toolboxes)
    if isOK(k)
        fprintf('%-14s OK\n', toolboxes{k});
    else
        fprintf('%-14s FAILED (see message above)\n', toolboxes{k});
    end
end

end

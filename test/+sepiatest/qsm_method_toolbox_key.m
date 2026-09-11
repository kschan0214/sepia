%% key = sepiatest.qsm_method_toolbox_key(qsmMethod)
%
% Description: maps a QSM dipole-inversion method name (as it appears in
% methodQSMName, sepia_universal_variables.m) to the field name in
% sepiatest.discover_toolboxes()'s output struct that gates whether it
% can run on this machine. Returns 'none' for methods that need no
% external toolbox, or 'unavailable' for methods this test suite cannot
% exercise at all yet (e.g. deep-learning methods with no distributed
% model checkpoint files in this environment) - callers should treat
% 'unavailable' as an unconditional skip.
%
function key = qsm_method_toolbox_key(qsmMethod)

switch qsmMethod
    case {'TKD', 'Closed-form solution', 'iLSQR', 'NDI'}
        key = 'none';
    case {'STI suite iLSQR', 'Star-QSM'}
        key = 'STISuite';
    case 'FANSI'
        key = 'FANSI';
    case 'MEDI'
        key = 'MEDI';
    case 'LSQR+HEIDI'
        key = 'HEIDI';
    case 'MRI Suscep. Calc.'
        key = 'MRISC';
    case {'Chi-separation', 'LPCNN', 'QSMnet+', 'xQSM'}
        % deep-learning / ONNX-checkpoint-backed methods: no checkpoint
        % files are distributed with SEPIA or tracked in this repo, so
        % this suite has no way to exercise them yet on a generic machine.
        key = 'unavailable';
    otherwise
        % unrecognised (e.g. a newly-added method not yet categorised
        % here): conservatively require nothing, but this should be
        % updated to the method's actual dependency.
        key = 'none';
end

end

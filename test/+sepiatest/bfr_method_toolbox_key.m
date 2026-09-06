%% key = sepiatest.bfr_method_toolbox_key(bfrMethod)
%
% Description: maps a background-field-removal method name (as it appears
% in methodBFRName, sepia_universal_variables.m) to the field name in
% sepiatest.discover_toolboxes()'s output struct that gates whether it
% can run on this machine. See sepiatest.qsm_method_toolbox_key.m for the
% 'none'/'unavailable' convention.
%
function key = bfr_method_toolbox_key(bfrMethod)

switch bfrMethod
    case 'VSHARP'
        key = 'none';
    case {'LBV', 'PDF', 'RESHARP', 'SHARP'}
        key = 'MEDI';
    case {'VSHARP (STI suite)', 'VSHARP (STI suite 2D)', 'iHARPERELLA'}
        key = 'STISuite';
    case 'BFRnet'
        % deep-learning method: no model checkpoint files distributed
        % with SEPIA or tracked in this repo
        key = 'unavailable';
    otherwise
        key = 'none';
end

end

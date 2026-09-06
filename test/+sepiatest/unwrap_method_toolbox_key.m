%% key = sepiatest.unwrap_method_toolbox_key(unwrapMethod)
%
% Description: maps a phase-unwrapping method name (as it appears in
% methodUnwrapName, sepia_universal_variables.m) to the field name in
% sepiatest.discover_toolboxes()'s output struct that gates whether it
% can run on this machine. See sepiatest.qsm_method_toolbox_key.m for the
% 'none'/'unavailable' convention; '3D best path' uses 'excluded' (a
% third, deliberate category - see test/README.md "Known exclusions"):
% it fills out-of-mask voxels with unseeded rand() and its own docstring
% says it only works on the DCCN cluster, so it is never run by this
% regression suite regardless of toolbox availability.
%
function key = unwrap_method_toolbox_key(unwrapMethod)

switch unwrapMethod
    case 'None'
        key = 'none';
    case {'Laplacian (MEDI)', 'Region growing (MEDI)', 'Graphcut'}
        key = 'MEDI';
    case 'Laplacian (STI suite)'
        key = 'STISuite';
    case 'ROMEO'
        key = 'MRITOOLS';
    case 'SEGUE'
        key = 'SEGUE';
    case '3D best path'
        key = 'excluded';
    otherwise
        key = 'none';
end

end

%% sepiatest.assume_toolbox_available(testCase, toolboxes, key, label)
%
% Description: shared assumeTrue/assumeFail gate used by every Tier-2
% matrix test before it runs a method. `key` is one of:
%   'none'        - method needs no external toolbox, always proceeds
%   'unavailable' - method has no model/checkpoint files in this
%                   environment (deep-learning add-ons); always skipped
%   'excluded'    - method is deliberately excluded from this regression
%                   suite (e.g. '3D best path' unwrap - nondeterministic
%                   and cluster-only, see test/README.md); always skipped
%   otherwise     - a field name in `toolboxes` (from
%                   sepiatest.discover_toolboxes()); skipped unless true
%
% Input
% --------------
% testCase  : a matlab.unittest.qualifications.Qualifiable
% toolboxes : struct from sepiatest.discover_toolboxes()
% key       : as returned by sepiatest.{qsm,bfr,unwrap}_method_toolbox_key
% label     : the method name, for the skip message
%
function assume_toolbox_available(testCase, toolboxes, key, label)

switch key
    case 'none'
        % always proceeds
    case 'unavailable'
        testCase.assumeFail(sprintf( ...
            '"%s" has no model/checkpoint files available in this environment - skipping.', label));
    case 'excluded'
        testCase.assumeFail(sprintf( ...
            '"%s" is deliberately excluded from this regression suite (see test/README.md) - skipping.', label));
    otherwise
        testCase.assumeTrue(isfield(toolboxes,key) && toolboxes.(key), ...
            sprintf('"%s" requires the %s toolbox, which is not installed on this machine - skipping.', label, key));
end

end

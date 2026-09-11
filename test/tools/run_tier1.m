function run_tier1()
%% run_tier1()
%
% Description: runs the Tier-1 toolbox-free smoke test suite
% (test/tier1_smoke/) and errors (non-zero exit under `matlab -batch`) if
% any test fails. Used both locally and by
% .github/workflows/tier1-smoke.yml.
%
import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.TAPPlugin
import matlab.unittest.plugins.ToFile

testRoot = fileparts(fileparts(mfilename('fullpath'))); % .../test
addpath(testRoot); % so the +sepiatest package resolves

suite  = TestSuite.fromFolder(fullfile(testRoot, 'tier1_smoke'), 'IncludingSubfolders', true);
runner = TestRunner.withTextOutput();

tapFile = fullfile(testRoot, 'tools', 'tier1-results.tap');
runner.addPlugin(TAPPlugin.producingVersion13(ToFile(tapFile)));

results = runner.run(suite);
disp(table(results));

% "Incomplete" (assumeTrue-filtered, e.g. no reference saved yet) is not a
% failure - only an actual Failed result should break the build.
assert(~any([results.Failed]), ...
    'run_tier1:testsFailed', 'One or more Tier-1 smoke tests failed - see output above.');

end

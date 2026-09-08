function run_tier2()
%% run_tier2()
%
% Description: runs the Tier-2 toolbox-dependent regression matrix
% (test/tier2_matrix/) against a small, deterministically-generated
% synthetic phantom (test/phantom/generate_synthetic_phantom.m) - unlike
% Tier 3, no real dataset is required. Rows requiring an unavailable
% toolbox, or missing a reference, are skipped (not failed). Errors
% (non-zero exit under `matlab -batch`) only if a row that actually ran
% failed its comparison.
%
import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.TAPPlugin
import matlab.unittest.plugins.ToFile

testRoot = fileparts(fileparts(mfilename('fullpath'))); % .../test
addpath(testRoot); % so the +sepiatest package resolves before class parsing

suite  = TestSuite.fromFolder(fullfile(testRoot, 'tier2_matrix'), 'IncludingSubfolders', true);
runner = TestRunner.withTextOutput();

tapFile = fullfile(testRoot, 'tools', 'tier2-results.tap');
runner.addPlugin(TAPPlugin.producingVersion13(ToFile(tapFile)));

results = runner.run(suite);
disp(table(results));

assert(~any([results.Failed]), ...
    'run_tier2:testsFailed', 'One or more Tier-2 regression tests failed - see output above.');

end

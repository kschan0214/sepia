%% p = sepiatest.test_root()
%
% Description: returns the absolute path of the test/ folder (parent of
% the +sepiatest package), used to locate test/references, test/tier1_smoke, etc.
%
function p = test_root()

thisFile   = mfilename('fullpath');   % .../test/+sepiatest/test_root
packageDir = fileparts(thisFile);     % .../test/+sepiatest
p          = fileparts(packageDir);   % .../test

end

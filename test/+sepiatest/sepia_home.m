%% p = sepiatest.sepia_home()
%
% Description: returns the absolute path of the SEPIA repository root
% (the directory containing sepia.m/sepiaIO.m), derived from this file's
% own location (test/+sepiatest/sepia_home.m is two levels below repo
% root) rather than assuming any particular current directory.
%
function p = sepia_home()

thisFile = mfilename('fullpath');           % .../test/+sepiatest/sepia_home
packageDir = fileparts(thisFile);           % .../test/+sepiatest
testDir    = fileparts(packageDir);         % .../test
p          = fileparts(testDir);            % repo root

end

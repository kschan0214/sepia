%% [RDF] = Wrapper_BFR_LBV(totalField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
%
% Input
% --------------
% totalField    : total field map (background + tissue fields), in Hz
% mask          : signal mask
% matrixSize    : size of the input image
% voxelSize     : spatial resolution of each dimension of the data, in mm
% algorParam    : structure contains fields with algorithm-specific parameter(s)
% headerAndExtraData : structure contains extra header info/data for the algorithm
%
% Output
% --------------
% RDF           : local field map
%
% Description: This is a wrapper function to access LBV for SEPIA
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 8 March 2020 (v0.8.0)
% Date last modified:
%
%
function [RDF] = Wrapper_BFR_LBV(totalField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
sepia_universal_variables;

% get algorithm parameters
algorParam  = check_and_set_algorithm_default(algorParam);
method      = algorParam.bfr.method;
tol         = algorParam.bfr.tol;
depth       = algorParam.bfr.depth;
peel        = algorParam.bfr.peel;

% add path
sepia_addpath('MEDI');

%% Display algorithm parameters
disp('The following parameters are being used...');
disp(['Tolerance    = ' num2str(tol)]);
disp(['Depth        = ' num2str(depth)]);
disp(['Peel         = ' num2str(peel)]);

%% main
RDF         = LBV(totalField,mask,matrixSize,voxelSize,tol,depth,peel);
deleteme    = dir('mask*.bin');
delete(deleteme(1).name);
%         system(['rm ' deleteme.folder filesep deleteme.name]);
       
end

%% set default parameter if not specified
function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;

try algorParam2.bfr.tol     = algorParam.bfr.tol;   catch; algorParam2.bfr.tol  = 1e-4; end
try algorParam2.bfr.depth   = algorParam.bfr.depth; catch; algorParam2.bfr.depth= 4;    end
try algorParam2.bfr.peel    = algorParam.bfr.peel;  catch; algorParam2.bfr.peel = 1;    end

end
%% function [chi, lambdaOptimal] = qsmIterativeLSQR(localField,mask,matrixSize,voxelSize,varargin)
%
% Description: compute QSM based on iterative LSQR
%
% Input
% _____
%   localField      : local field perturbatios
%   mask            : user-defined mask
%   matrixSize      : image matrix size
%   voxelSize       : spatial resolution of image 
%   varargin        : flags with
%       'lambda'    -   user define regularisation parameter
%       'tol'       -	error tolerance
%       'iteration' -   no. of maximum iteration for iLSQR
%       'weight'    -   weighting of error computation
%       'initGuess' -   initial guess for iLSQR
%       'optimise'  -	self-define regularisation based on curvature of 
%                       L-curve
% 
% Ouput
% _____
%   chi             : QSM
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 24 March 2017
% Date last modified: 6 September 2017
%
function [chi, lambdaOptimal] = qsmIterativeLSQR(localField,mask,matrixSize,voxelSize,varargin)
lambdaOptimal = [];
%% Parsing varargin
[lambda, tol, maxiter, wmap, initGuess, optimise] = parse_varargin_iLSQR(varargin);

% display message
disp('The following parameters are used:');
fprintf('Regularisation parameter: %f \n',lambda);
fprintf('Maximum iterations: %i \n',maxiter);
fprintf('Tolerance: %f \n',tol);

if isempty(initGuess)
    initGuess = zeros(matrixSize);
    disp('Initial guesses = 0 ');
end
if isempty(wmap)
    wmap = ones(matrixSize);
    disp('No weighting');
end

% dipole kernel
kernel = DipoleKernel(matrixSize,voxelSize);

%% Core
if optimise
    [~, lambdaOptimal] = qsmClosedFormL2(localField,mask,matrixSize,voxelSize,'optimise',optimise);
    lambda = lambdaOptimal;
end
% defining gradient operators in k-space
[Ex,Ey,Ez,~] = GradientOperatorKspace(matrixSize);

% input parameters for lsqr optimization
% optimization parameters
params_in.msk = wmap;
params_in.Ex = Ex;
params_in.Ey = Ey;
params_in.Ez = Ez;
params_in.kernel = kernel;
params_in.dims = matrixSize;
params_in.lambda = lambda;

% right hand of the system of linear equations  Ax=b 
% KC: b is the right side 4 submatrices in r-space
b = cat( 4 , localField.*params_in.msk , zeros(matrixSize) , zeros(matrixSize) , zeros(matrixSize) );
b = b(:);

% initial guess of susceptibility distribution
% x = zeros( prod(params_in.dims) , 1);

% using lsqr to obtain the solution
tic
[chi_L2iterative,~,~,~] = lsqr( @(x,tflag) QSM_gradRegul(x,params_in,tflag),b,tol,maxiter,[],[],initGuess(:));
toc

% reshaping into a 3D
chi = reshape( chi_L2iterative , params_in.dims ).*mask;
end

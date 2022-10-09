%% function [chi, lambdaOptimal] = qsmIterativeLSQR(localField,mask,matrixSize,voxelSize,varargin)
%
% Description: compute QSM based on iterative LSQR
%
% Input
% _____
%   localField      : local field perturbations
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
%       'b0dir'     -   B0 direction
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
[lambda,tol,maxiter,wmap,initGuess,optimise,b0dir] = parse_varargin_iLSQR(varargin);

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
kernel = DipoleKernel(matrixSize,voxelSize,b0dir);

%% Core
if optimise
    [~, lambdaOptimal] = qsmClosedFormL2(localField,mask,matrixSize,voxelSize,'optimise',optimise,'b0dir',b0dir);
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

%% parser
function [lambda,tol,maxiter,wmap,initGuess,optimise,b0dir] = parse_varargin_iLSQR(arg)

% predefine parameters
lambda      = 1e-1;
tol         = 1e-3;
maxiter     = 50;
wmap        = [];
initGuess   = [];
optimise    = false;
b0dir       = [0,0,1];

% use user defined input if any
if ~isempty(arg)
    for kvar = 1:length(arg)
        if strcmpi(arg{kvar},'lambda')
            lambda = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'tol')
            tol = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'iteration')
            maxiter = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'weight')
            wmap = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'initGuess')
            initGuess = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'optimise')
            optimise = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'b0dir')
            b0dir = arg{kvar+1};
        end
    end
end

end

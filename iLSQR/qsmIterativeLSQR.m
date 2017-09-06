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
% Date last modified: 28 June 2017
%
function [chi, lambdaOptimal] = qsmIterativeLSQR(localField,mask,matrixSize,voxelSize,varargin)
lambdaOptimal = [];
%% Parsing varargin
[lambda, tol, maxiter, wmap, initGuess, optimise] = parse_vararginiLSQR(varargin);

if isempty(initGuess)
    initGuess = zeros(matrixSize);
end
if isempty(wmap)
    wmap = ones(matrixSize);
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
% magn = imDataParameter.magn;
% params_in.msk = single(mask.*magn)./norm(single(mask(:).*magn(:)));
% params_in.msk = single(mask.*magn)./max(single(mask(:).*magn(:)));
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
[chi_L2iterative,FLAG,RELRES,ITER] = lsqr( @(x,tflag) QSM_gradRegul(x,params_in,tflag),b,tol,maxiter,[],[],initGuess(:));
toc

% reshaping into a 3D
chi = reshape( chi_L2iterative , params_in.dims ).*mask;
end

%% Parsing varargin
% function [lambda, tol, maxiter, wmap, initGuess, optimise] = parse_vararginiLSQR(arg)
% lambda = 1e-1;
% tol = 1e-3;
% maxiter = 50;
% wmap = [];
% initGuess = [];
% optimise = false;
% 
% if ~isempty(arg)
%     for kvar = 1:length(arg)
%         if strcmpi(arg{kvar},'lambda')
%             lambda = arg{kvar+1};
%         end
%         if strcmpi(arg{kvar},'tol')
%             tol = arg{kvar+1};
%         end
%         if strcmpi(arg{kvar},'iteration')
%             maxiter = arg{kvar+1};
%         end
%         if strcmpi(arg{kvar},'weight')
%             wmap = arg{kvar+1};
%         end
%         if strcmpi(arg{kvar},'initGuess')
%             initGuess = arg{kvar+1};
%         end
%         if strcmpi(arg{kvar},'optimise')
%             optimise = arg{kvar+1};
%         end
%     end
% end
% end
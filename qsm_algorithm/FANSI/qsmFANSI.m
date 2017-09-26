%% function chi = qsmFANSI(localField,mask,matrixSize,voxelSize,varargin)
%
% Usage:
%   chi = qsmFANSI(localField,mask,matrixSize,voxelSize,...
%           'tol',1,'lambda',3e-5,'mu',5e-5,'iteration',50,'weight',wmap)
%   chi = qsmFANSI(localField,mask,matrixSize,voxelSize,...
%           'tol',1,'lambda',3e-5,'mu',5e-5,'iteration',50,'weight',wmap,...
%           'linear','tv');
%
% Input
% --------------
%   localField      : local field perturbations
%   mask            : user-defined mask
%   matrixSize      : image matrix size
%   voxelSize       : spatial resolution of image 
%   varargin        : flags with
%       'lambda'    -   user defined regularisation parameter for gradient L1 penalty
%       'mu'        -	user defined regularisation parameter for gradient consistency 
%       'tol'       - tolerance for iteration
%       'iteration' - maximum number of iterations
%       'weight'    - weighting of error computation, values have to be <1
%       'linear'    - linear solver
%       'nonlinear' - nonlinear solver
%       'tv'        - Total variation constraints
%       'tgv'       - Total general variation constraints
%       'b0dir'     - B0 direction
%
% Output
% --------------
%   chi             : QSM
%
% Description: QSM computation based on linear/nonlinear TV/TGV constraints
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 14 July 2017
% Date last modified: 26 September 2017
%
%
function chi = qsmFANSI(localField,mask,matrixSize,voxelSize,varargin)
% parse input argument
[mu1,alpha1,tol,maxiter,wmap,solver,constraint,b0dir] = parse_varargin_FANSI(varargin);

% display message
fprintf('Regularisation parameter for gradient L1: %f \n',alpha1);
fprintf('Regularisation parameter for gradient consistency: %f \n',mu1);
fprintf('Tolerance: %f \n',tol);
fprintf('Maximum iterations: %i \n',maxiter);

if isempty(wmap)
    wmap = ones(matrixSize);
end
params = [];
params.K = DipoleKernel(matrixSize,voxelSize,b0dir);
params.N = matrixSize;
params.input = localField;
params.weight = wmap; 

params.maxOuterIter = maxiter;

params.mu1 = mu1;          % gradient consistency
params.alpha1 = alpha1;       % gradient L1 penalty (TGV constraint)
params.tol_update = tol;

switch solver
    case 'linear'
        disp('Solver: Linear');
        switch constraint
            case 'TV'
                disp('Constraint: TV');
                out = wTV(params);
            case 'TGV'
                disp('Constraint: TGV');
                out = wTGV(params);
        end
    case 'nonlinear'
        disp('Solver: Non-linear');
        switch constraint
            case 'TV'
                disp('Constraint: TV');
                out = nlTV(params);
            case 'TGV'
                disp('Constraint: TGV');
                out = nlTGV(params);
        end
end

chi = real(out.x).*mask;

end

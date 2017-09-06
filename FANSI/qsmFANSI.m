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
%       'weight'    - weighting of error computation
%       'linear'    - linear solver
%       'nonlinear' - nonlinear solver
%       'tv'        - Total variation constraints
%       'tgv'        - Total general variation constraints
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
% Date last modified: 6 September 2017
%
%
function chi = qsmFANSI(localField,mask,matrixSize,voxelSize,varargin)
% parse input argument
[mu1,alpha1,tol,maxiter,wmap,solver,constraint] = parse_varargin_FANSI(varargin);
params = [];
params.K = DipoleKernel(matrixSize,voxelSize);
params.N = matrixSize;
params.input = localField;
params.weight = wmap; 

params.maxOuterIter = maxiter;

params.mu1 = mu1;          % gradient consistency
params.alpha1 = alpha1;       % gradient L1 penalty (TGV constraint)
params.tol_update = tol;

switch solver
    case 'linear'
        switch constraint
            case 'TV'
                out = wTV(params);
            case 'TGV'
                out = wTGV(params);
        end
    case 'nonlinear'
        switch constraint
            case 'TV'
                out = nlTV(params);
            case 'TGV'
                out = nlTGV(params);
        end
end

chi = real(out.x).*mask;

end

%% parse argument input
% function [mu1,alpha1,tol,maxiter,wmap,solver,constraint]=parse_vararginFANSI(arg)
% alpha1 = 3e-5;
% mu1 = 5e-5;
% maxiter = 40;
% wmap = [];
% solver = 'nonlinear';
% constraint = 'TGV';
% tol = 1;
% 
% if ~isempty(arg)
%     for kvar = 1:length(arg)
%         if strcmpi(arg{kvar},'lambda')
%             alpha1 = arg{kvar+1};
%         end
%         if strcmpi(arg{kvar},'mu')
%             mu1 = arg{kvar+1};
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
%         if strcmpi(arg{kvar},'linear')
%             solver = 'linear';
%         end
%         if strcmpi(arg{kvar},'nonlinear')
%             solver = 'nonlinear';
%         end
%         if strcmpi(arg{kvar},'tv')
%             constraint = 'TV';
%         end
%         if strcmpi(arg{kvar},'tgv')
%             constraint = 'TGV';
%         end
%     end
% end
% end
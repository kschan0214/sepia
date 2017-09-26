%% function [chi, lamdaOptimal] = qsmMacro(localField,mask,matrixSize,voxelSize,varargin)
%
% Usage: 
%       chi = qsmMacro(localField,mask,matrixSize,voxelSize,...
%               'method','TKD','threshold',0.15);
%       chi = qsmMacro(localField,mask,matrixSize,voxelSize,...
%               'method','ClosedFormL2','lambda',0.1,'optimise',false);
%       chi = qsmMacro(localField,mask,matrixSize,voxelSize,...
%               'method','iLSQR','lambda',0.1,'optimise',false,'tol',1e-3,...
%               'iteration',100,'weight',wmap,'initGuess',initGuessmap);
%       chi = qsmMacro(localField,mask,matrixSize,voxelSize,...
%               'method','STISuiteiLSQR','threshold',0.01,'iteration',100,...
%               'tol_step1',0.01,'tol_step2',0.001,'b0dir',[0,0,1],'TE',1,...
%               'B0',3,'padsize',[4,4,4]);
%       chi = qsmMacro(localField,mask,matrixSize,voxelSize,...
%               'method','FANSI','tol',1,'lambda',3e-5,'mu',5e-5,'iteration',50,'weight',wmap,...
%               'tgv','nonlinear');
%       chi = qsmMacro(totalField,mask,matrixSize,voxelSize,...
%               'method','SSVSHARP','tol',tol,'lambda',lambda,'iteration',maxiter,'magnitude',magn,...
%               'vkernel',Kernel_Sizes);
%
% Description: Wrapper for QSM inversion problem (default using TKD)
%   Flags:
%       'method'        : QSM method, 
%                          'TKD', 'ClosedFormL2', 'iLSQR', 'FANSI',
%                          'ssvsharp', 'STISuiteiLSQR'
%       'b0dir'         : B0 direction
%
%       TKD
%       ----------------
%       'threshold'     : threshold for TKD
%
%       ClosedFormL2
%       ----------------           
%       'lambda'        : regularisation value
%       'optimise'      : self-define regularisation based on curvature of 
%                         L-curve
%
%       iLSQR
%       ----------------
%       'lambda'        : regularisation value
%       'tol'           : error tolerance
%       'iteration' 	: no. of maximum iteration for iLSQR
%       'weight'    	: weighting of error computation
%       'initGuess'     : initial guess for iLSQR
%       'optimise'      : self-define regularisation based on curvature of 
%                         L-curve
%
%       STISuiteiLSQR
%       ----------------
%       'threshold'     : threshold for STISuiteiLSQR
%       'iteration' 	: no. of maximum iteration for iLSQR
%       'tol_step1'     : error tolerance
%       'tol_step2'     : error tolerance
%       'b0dir'         : main magnetic field direction
%       'TE'            : echo time
%       'fieldStrength' : magntic field strength of the scanner
%       'padsize'       : size for padarray to increase numerical accuracy
%
%       FANSI
%       ----------------
%       'lambda'        : user defined regularisation parameter for gradient L1 penalty
%       'mu'            : user defined regularisation parameter for gradient consistency 
%       'tol'           : tolerance for iteration
%       'iteration'     : maximum number of iterations
%       'weight'        : weighting of error computation
%       'linear'        : linear solver
%       'nonlinear'     : nonlinear solver
%       'tgv'           : Total generalised variation constraint
%       'tv'            : Total variation constraints
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 28 June 2017
% Date last modified: 6 September 2017
%
function [chi, lamdaOptimal] = qsmMacro(localField,mask,matrixSize,voxelSize,varargin)
lamdaOptimal = [];
%% Parsing argument input flags
if ~isempty(varargin)
    for kvar = 1:length(varargin)
        if strcmpi(varargin{kvar},'method')
            switch lower(varargin{kvar+1})
                case 'tkd'
                    method = 'TKD';
                    [thre_tkd,b0dir] = parse_varargin_TKD(varargin);
                    break
                case 'closedforml2'
                    method = 'CFL2';
                    [lambda,optimise,b0dir] = parse_varargin_CFL2norm(varargin);
                    break
                case 'ilsqr'
                    method = 'iLSQR';
                    [lambda, tol, maxiter, wmap, initGuess, optimise,b0dir] = parse_varargin_iLSQR(varargin);
                    if isempty(wmap)
                        wmap = ones(matrixSize);
                    end
                    if isempty(initGuess)
                        initGuess = zeros(matrixSize);
                    end
                    break
                case 'stisuiteilsqr'
                    method = 'STISuiteiLSQR';
                    algoPara = parse_varargin_STISuiteiLSQR(varargin);
                    algoPara.voxelsize= voxelSize(:).';
                    break
                case 'fansi'
                    method = 'FANSI';
                    [mu1,alpha1,tol,maxiter,wmap,solver,constraint,b0dir]=parse_varargin_FANSI(varargin);
                case 'ssvsharp'
                    method = 'SSVSHARP';
                    [lambda,magn,tol,maxiter,Kernel_Sizes,b0dir]=parse_varargin_SSQSM(varargin);
            end
        end
    end
else
    % predefine paramater: if no varargin, use TKD
    disp('No method selected. Using default setting...');
    method = 'TKD';
    thre_tkd = 0.15;
end

disp(['The following QSM algorithm will be used: ' method]);

%% qsm algorithm
switch method
    case 'TKD'
        chi = qsmTKD(localField,mask,matrixSize,voxelSize,'threshold',thre_tkd,'b0dir',b0dir);
    case 'CFL2'
        [chi, lamdaOptimal] = qsmClosedFormL2(localField,mask,matrixSize,voxelSize,...
            'lambda',lambda,'optimise',optimise,'b0dir',b0dir);
    case 'iLSQR'
        chi = qsmIterativeLSQR(localField,mask,matrixSize,voxelSize,...
            'lambda',lambda,'tol',tol,'iteration',maxiter,'weight',wmap,...
            'initGuess',initGuess,'optimise',optimise,'b0dir',b0dir);
    case 'STISuiteiLSQR'
        chi = QSM_iLSQR(localField,mask,'params',algoPara);
    case 'FANSI'
        chi = qsmFANSI(localField,mask,matrixSize,voxelSize,...
          'tol',tol,'lambda',alpha1,'mu',mu1,'iteration',maxiter,'weight',wmap,...
          solver,constraint,'b0dir',b0dir);
    case 'SSVSHARP'
        chi = qsmSingleStepVSHARP(localField,mask,matrixSize,voxelSize,...
            'tol',tol,'lambda',lambda,'iteration',maxiter,'magnitude',magn,...
            'b0dir',b0dir,'vkernel',Kernel_Sizes);
end

end

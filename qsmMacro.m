%% function [chi, lamdaOptimal] = qsmMacro(localField,mask,matrixSize,voxelSize,varargin)
%
% Usage: 
%       chi = qsmMacro(localField,mask,matrixSize,voxelSize,...
%               'method','TKD','threshold',0.15);
%       chi = qsmMacro(localField,mask,matrixSize,voxelSize,...
%               'method','ClosedFormL2','lambda',0.1,'optimise',false);
%       chi = qsmMacro(localField,mask,matrixSize,voxelSize,...
%               'method','iLSQR','lambda',0.1,'optimise',false,'tol',1e-3,..
%               'iteration',100,'weight',wmap,'initGuess',initGuessmap);
%       chi = qsmMacro(localField,mask,matrixSize,voxelSize,...
%               'method','STISuiteiLSQR','threshold',0.01,'iteration',100,
%               'tol_step1',0.01,'tol_step2',0.001,'b0dir',[0,0,1],'TE',1
%               'fieldStrength',3,'padsize',[4,4,4]);
%
% Description: Wrapper for QSM (default using TKD)
%   Flags:
%       'method'        : QSM method, 
%                          'TKD', 'ClosedFormL2', 'iLSQR' and
%                          'STISuiteiLSQR'
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
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 28 June 2017
% Date last modified:
%
function [chi, lamdaOptimal] = qsmMacro(localField,mask,matrixSize,voxelSize,varargin)
lamdaOptimal = [];
%% Parsing argument input flags
if ~isempty(varargin)
    for kvar = 1:length(varargin)
        if strcmpi(varargin{kvar},'method')
            switch lower(varargin{kvar+1})
                case 'tkd'
                    [method, thre_tkd] = parse_vararginTKD(varargin);
                    break
                case 'closedforml2'
                    [method, lambda, optimise] = parse_vararginCFL2(varargin);
                    break
                case 'ilsqr'
                    [method, lambda, tol, maxiter, wmap, initGuess, optimise] = parse_vararginiLSQR(varargin);
                    if isempty(wmap)
                        wmap = ones(matrixSize);
                    end
                    if isempty(initGuess)
                        initGuess = zeros(matrixSize);
                    end
                    break
                case 'stisuiteilsqr'
                    [method, algoPara] = parse_vararginSTISuiteiLSQR(varargin);
                    algoPara.voxelsize= voxelSize;
                    break
            end
        end
    end
else
    % predefine paramater: if no varargin, use LBV
    disp('No method selected. Using the default setting:');
    method = 'TKD'
    thre_tkd = 0.15
end

%% qsm algorithm
switch method
    case 'TKD'
        chi = qsmTKD(localField,mask,matrixSize,voxelSize,'threshold',thre_tkd);
    case 'CFL2'
        [chi, lamdaOptimal] = qsmClosedFormL2(localField,mask,matrixSize,voxelSize,...
            'lambda',lambda,'optimise',optimise);
    case 'iLSQR'
        chi = qsmIterativeLSQR(localField,mask,matrixSize,voxelSize,...
            'lambda',lambda,'tol',tol,'iteration',maxiter,'weight',wmap,...
            'initGuess',initGuess,'optimise',optimise);
    case 'STISuiteiLSQR'
        chi = QSM_iLSQR(localField,mask,'params',algoPara);
end

end

%% Parsing varargin
% TKD
function [method, thre_tkd] = parse_vararginTKD(arg)
method = 'TKD';
thre_tkd = 0.15;
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'threshold')
        thre_tkd = arg{kkvar+1};
        continue
    end
end
end

% Closed form solution L2norm
function [method, lambda, optimise] = parse_vararginCFL2(arg)
method = 'CFL2';
lambda = 1e-1;
optimise = false;
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'lambda')
        lambda = arg{kkvar+1};
        continue
    end
    if  strcmpi(arg{kkvar},'optimise')
        optimise = arg{kkvar+1};
        continue
    end
end
end

% iLSQR
function [method, lambda, tol, maxiter, wmap, initGuess, optimise] = parse_vararginiLSQR(arg)
method = 'iLSQR';
lambda = 1e-1;
tol = 1e-3;
maxiter = 50;
wmap = [];
initGuess = [];
optimise = false;
for kvar = 1:length(arg)
    if strcmpi(arg{kvar},'lambda')
        lambda = arg{kvar+1};
        continue
    end
    if strcmpi(arg{kvar},'tol')
        tol = arg{kvar+1};
        continue
    end
    if strcmpi(arg{kvar},'iteration')
        maxiter = arg{kvar+1};
        continue
    end
    if strcmpi(arg{kvar},'weight')
        wmap = arg{kvar+1};
        continue
    end
    if strcmpi(arg{kvar},'initGuess')
        initGuess = arg{kvar+1};
        continue
    end
    if  strcmpi(arg{kkvar},'optimise')
        optimise = arg{kkvar+1};
        continue
    end
end
end

% STI suite iLSQR
function [method, params] = parse_vararginSTISuiteiLSQR(arg)
method = 'STISuiteiLSQR';
params.H=[0 0 1];
params.niter=100;
params.TE=1;
params.B0=3;
params.tol_step1=0.01;
params.tol_step2=0.001;
params.Kthreshold=0.25;
params.padsize=[4,4,4];

for kvar = 1:length(arg)
    if strcmpi(arg{kvar},'b0dir')
        params.H = arg{kvar+1};
        continue
    end
    if strcmpi(arg{kvar},'tol_step1')
        params.tol_step1 = arg{kvar+1};
        continue
    end
    if strcmpi(arg{kvar},'tol_step2')
        params.tol_step2 = arg{kvar+1};
        continue
    end
    if strcmpi(arg{kvar},'iteration')
        params.niter = arg{kvar+1};
        continue
    end
    if strcmpi(arg{kvar},'TE')
        params.TE = arg{kvar+1};
        continue
    end
    if strcmpi(arg{kvar},'fieldStrength')
        params.B0 = arg{kvar+1};
        continue
    end
    if strcmpi(arg{kvar},'threshold')
        params.Kthreshold = arg{kvar+1};
        continue
    end
    if strcmpi(arg{kvar},'padsize')
        params.padsize = arg{kvar+1};
        continue
    end
end
end
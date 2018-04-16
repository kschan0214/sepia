%% function output = function_name(input)
%
% Usage:
%
% Input
% --------------
%
% Output
% --------------
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 
% Date last modified:
%
%
function [isBET,mask,unwrap,subsampling,BFR,refine,BFR_tol,BFR_depth,BFR_peel,BFR_iteration,...
    BFR_padSize,BFR_radius,BFR_alpha,BFR_threshold,QSM_method,QSM_threshold,QSM_lambda,...
    QSM_optimise,QSM_tol,QSM_maxiter,QSM_tol1,QSM_tol2,QSM_padsize,QSM_mu1,QSM_solver,QSM_constraint,...
    exclude_threshold,QSM_radius,QSM_zeropad,QSM_wData,QSM_wGradient,QSM_lambdaCSF,QSM_isSMV,QSM_merit,isEddyCorrect] = parse_varargin_QSMHub(arg)

mask=[];
isBET=false;
unwrap='Laplacian';subsampling=1;
BFR='LBV';refine=false;BFR_tol=1e-4;BFR_depth=4;BFR_peel=2;BFR_iteration=50;
BFR_padSize=40;BFR_radius=4;BFR_alpha=0.01;BFR_threshold=0.03;
QSM_method='TKD';QSM_threshold=0.15;QSM_lambda=0.13;QSM_optimise=false;
QSM_tol=1e-3;QSM_maxiter=50;QSM_tol1=0.01;QSM_tol2=0.001;QSM_padsize=[4,4,4];
QSM_mu1=5e-5;QSM_solver='linear';QSM_constraint='tv';
exclude_threshold = 1;
QSM_radius=5;QSM_zeropad=0;QSM_wData=1;QSM_wGradient=1;QSM_lambdaCSF=100;
QSM_isSMV=false;QSM_merit=false;isEddyCorrect=false;

if ~isempty(arg)
    for kvar = 1:length(arg)
        if strcmpi(arg{kvar},'FSLBet')
            isBET = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'mask')
            mask = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'unwrap')
            unwrap = arg{kvar+1};
        end
%         if strcmpi(arg{kvar},'unit')
%             unit = arg{kvar+1};
%         end
        if  strcmpi(arg{kvar},'Subsampling')
            subsampling = arg{kvar+1};
            continue
        end
        if  strcmpi(arg{kvar},'BFR')
            BFR = arg{kvar+1};
            continue
        end
        if  strcmpi(arg{kvar},'refine')
            refine = arg{kvar+1};
            continue
        end
        if strcmpi(arg{kvar},'BFR_tol')
            BFR_tol = arg{kvar+1};
            continue
        end
        if  strcmpi(arg{kvar},'depth')
            BFR_depth = arg{kvar+1};
            continue
        end
        if strcmpi(arg{kvar},'peel')
            BFR_peel = arg{kvar+1};
            continue
        end
        if strcmpi(arg{kvar},'BFR_iteration')
            BFR_iteration = arg{kvar+1};
            continue
        end
        if strcmpi(arg{kvar},'BFR_padsize')
            BFR_padSize = arg{kvar+1};
            continue
        end
        if strcmpi(arg{kvar},'BFR_radius')
            BFR_radius = arg{kvar+1};
            continue
        end
        if  strcmpi(arg{kvar},'BFR_alpha')
            BFR_alpha = arg{kvar+1};
            continue
        end
        if  strcmpi(arg{kvar},'BFR_threshold')
            BFR_threshold = arg{kvar+1};
            continue
        end
        if  strcmpi(arg{kvar},'QSM')
            QSM_method = arg{kvar+1};
            continue
        end
        if strcmpi(arg{kvar},'QSM_threshold')
            QSM_threshold = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'QSM_lambda')
            QSM_lambda = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'QSM_optimise')
            QSM_optimise = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'QSM_tol')
            QSM_tol = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'QSM_iteration')
            QSM_maxiter = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'QSM_tol1')
            QSM_tol1 = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'QSM_tol2')
            QSM_tol2 = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'QSM_padsize')
            QSM_padsize = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'QSM_mu')
            QSM_mu1 = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'linear')
            QSM_solver = 'linear';
        end
        if strcmpi(arg{kvar},'non-linear')
            QSM_solver = 'nonlinear';
        end
        if strcmpi(arg{kvar},'tv')
            QSM_constraint = 'TV';
        end
        if strcmpi(arg{kvar},'tgv')
            QSM_constraint = 'TGV';
        end
        if strcmpi(arg{kvar},'exclude_threshold')
            exclude_threshold = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'QSM_zeropad')
            QSM_zeropad = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'QSM_wData')
            QSM_wData = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'QSM_wGradient')
            QSM_wGradient = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'QSM_radius')
            QSM_radius = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'QSM_isSMV')
            QSM_isSMV = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'QSM_merit')
            QSM_merit = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'QSM_isLambdaCSF')
            QSM_isLambdaCSF = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'QSM_lambdaCSF')
            QSM_lambdaCSF = arg{kvar+1};
        end
       if strcmpi(arg{kvar},'eddy')
           isEddyCorrect = arg{kvar+1};
       end
%         if strcmpi(arg{kvar},'QSM_weight')
%             QSM_wmap = arg{kvar+1};
%         end
    end
end
end
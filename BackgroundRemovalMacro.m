%% RDF = BackgroundRemovalMacro(bkgField,mask,matrixSize,voxelSize,varargin)
%
% Usage: 
%       RDF = BackgroundRemovalMacro(fieldMap,mask,matrixSize,voxelSize,...
%               'method','LBV','refine',true,'tol',1e-4,'depth',4,'peel',2);
%       RDF = BackgroundRemovalMacro(fieldMap,mask,matrixSize,voxelSize,...
%               'method','RDF','refine',true,'tol',0.1,'b0dir',[0,0,1],...
%               'iteration',100,'CGsolver',true,'noisestd',N_std);
%       RDF = BackgroundRemovalMacro(fieldMap,mask,matrixSize,voxelSize,...
%               'method','SHARP','refine',true,'radius',4,'threshold',0.03);
%       RDF = BackgroundRemovalMacro(fieldMap,mask,matrixSize,voxelSize,...
%               'method','RESHARP','refine',true,'radius',4,'threshold',0.01);
%       RDF = BackgroundRemovalMacro(fieldMap,mask,matrixSize,voxelSize,...
%               'method','VSHARP','refine',true); 
%       RDF = BackgroundRemovalMacro(fieldMap,mask,matrixSize,voxelSize,...
%               'method','iHARPERELLA','refine',true,'iteration',100);
%
% Description: Wrapper for background field removal (default using LBV)
%   Flags:
%       'method'        : background revomal method, 
%                          'LBV', 'PDF', 'SHARP', 'RESHARP', 'VSHARP, and
%                          'iHARPERELLA'
%       'refine'        : refine the RDF by 5th order polynomial fitting
%
%       LBV
%       ----------------
%       'tol'           : error tolerance
%       'depth'         : multigrid level
%       'peel'          : thickness of the boundary layer to be peeled off
%
%       PDF
%       ----------------
%       'tol'           : error tolerance
%       'b0dir'         : B0 direction (e.g. [x,y,z])
%       'iteration'     : no. of maximum iteration for cgsolver
%       'CGsolver'      : true for default cgsolver; false for matlab pcg
%                         solver
%       'noisestd'      : noise standard deviation on the field map
%
%       SHARP
%       ----------------
%       'radius'        : radius of the spherical mean value operation
%       'threshold'     : threshold used in Truncated SVD (SHARP)
%
%       RESHARP
%       ----------------
%       'radius'        : radius of the spherical mean value operation
%       'alpha'         : regularizaiton parameter used in Tikhonov
%
%       VSHARP
%       ----------------
%
%       iHARPERELLA
%       ----------------
%       'iteration'     : no. of maximum iteration
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 28 June 2017
% Date last modified:
%
function RDF = BackgroundRemovalMacro(bkgField,mask,matrixSize,voxelSize,varargin)
%% Parsing argument input flags
if ~isempty(varargin)
    for kvar = 1:length(varargin)
        if strcmpi(varargin{kvar},'method')
            switch lower(varargin{kvar+1})
                case 'lbv'
                    [method, tol, depth, peel, refine] = parse_vararginLBV(varargin);
                    break
                case 'pdf'
                    [method, B0_dir, tol, iteration, CGdefault, N_std, refine] = parse_vararginPDF(matrixSize,varargin);
                    break
                case 'sharp'
                    [method, radius, threshold, refine] = parse_vararginSHARP(varargin);
                    break
                case 'resharp'
                    [method, radius, alpha, refine] = parse_vararginRESHARP(varargin);
                    break
                case 'vsharp'
                    [method, refine] = parse_vararginVSHARP(varargin);
                    break
                case 'iharperella'
                    [method, iteration, refine] = parse_vararginiHARPERELLA(varargin);
                    break
            end
        end
    end
else
    % predefine paramater: if no varargin, use LBV
    disp('No method selected. Using the default setting:');
    method = 'LBV'
    tol = 1e-4
    depth = 4
    peel = 1
    refine = false
end

%% background field removal
switch method
    case 'LBV'
        RDF = LBV(bkgField,mask,matrixSize,voxelSize,tol,depth,peel);
    case 'PDF'
        RDF = PDF(bkgField, N_std, mask,matrixSize,voxelSize, B0_dir, ...
            'tol', tol,'iteration', iteration,'CGsolver', CGdefault);
    case 'SHARP'
        RDF = SHARP(bkgField, mask, matrixSize, voxelSize, radius,threshold);
    case 'RESHARP'
        RDF = RESHARP(bkgField, mask, matrixSize, voxelSize, radius, alpha);
    case 'VSHARP'
        RDF = V_SHARP(bkgField, mask,'voxelsize',voxelSize);
    case 'iHARPERELLA'
        RDF = iHARPERELLA(bkgField, mask,'voxelsize',voxelSize,'niter',iteration);
end

%% If refine is needed, do it now
if refine
     [~,RDF,~]=PolyFit(RDF,RDF~=0,5);
end

end

%% Parsing varargin
% LBV
function [method, tol, depth, peel, refine] = parse_vararginLBV(arg)
method = 'LBV';
tol = 1e-4;
depth = 4;
peel = 1;
refine = false;
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'tol')
        tol = arg{kkvar+1};
        continue
    end
    if  strcmpi(arg{kkvar},'depth')
        depth = arg{kkvar+1};
        continue
    end
    if strcmpi(arg{kkvar},'peel')
        peel = arg{kkvar+1};
        continue
    end
    if strcmpi(arg{kkvar},'refine')
        refine = arg{kkvar+1};
        continue
    end
end
end

% PDF
function [method, B0_dir, tol, iteration, CGdefault, N_std, refine] = parse_vararginPDF(matrixSize,arg)
method = 'PDF';
B0_dir = [0,0,1];
tol = 0.1;
iteration = 30;
CGdefault = true;
N_std = ones(matrixSize)*1e-4;
refine = false;
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'b0dir')
        B0_dir = arg{kkvar+1};
        continue
    end
    if  strcmpi(arg{kkvar},'tol')
        tol = arg{kkvar+1};
        continue
    end
    if strcmpi(arg{kkvar},'iteration')
        iteration = arg{kkvar+1};
        continue
    end
    if strcmpi(arg{kkvar},'CGsolver')
        CGdefault = arg{kkvar+1};
        continue
    end
    if strcmpi(arg{kkvar},'noisestd')
        N_std = arg{kkvar+1};
        continue
    end
    if strcmpi(arg{kkvar},'refine')
        refine = arg{kkvar+1};
        continue
    end
end
end

% SHARP
function [method, radius, threshold, refine] = parse_vararginSHARP(arg)
method = 'SHARP';
radius = 4;
threshold = 0.03;
refine = false;
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'radius')
        radius = arg{kkvar+1};
        continue
    end
    if  strcmpi(arg{kkvar},'threshold')
        threshold = arg{kkvar+1};
        continue
    end
    if strcmpi(arg{kkvar},'refine')
        refine = arg{kkvar+1};
        continue
    end
end
end

% RESHARP
function [method, radius, alpha, refine] = parse_vararginRESHARP(arg)
method = 'RESHARP';
radius = 4;
alpha = 0.01;
refine = false;
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'radius')
        radius = arg{kkvar+1};
        continue
    end
    if  strcmpi(arg{kkvar},'alpha')
        alpha = arg{kkvar+1};
        continue
    end
    if strcmpi(arg{kkvar},'refine')
        refine = arg{kkvar+1};
        continue
    end
end
end

% VSHARP
function [method, refine] = parse_vararginVSHARP(arg)
method = 'VSHARP';
refine = false;
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'refine')
        refine = arg{kkvar+1};
        continue
    end
end
end

% iHARPERELLA
function [method, iteration, refine] = parse_vararginiHARPERELLA(arg)
method = 'iHARPERELLA';
iteration = 100;
refine = false;
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'iteration')
        iteration = arg{kkvar+1};
        continue
    end
    if strcmpi(arg{kkvar},'refine')
        refine = arg{kkvar+1};
        continue
    end
end
end
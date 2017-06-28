%% RDF = BackgroundRemovalMacro(bkgField,mask,matrixSize,voxelSize,varargin)
%
% Description: Wrapper for background field removal (default using LBV)
%   Flags:
%       'method'        : background revomal method, 
%                          'LBV', 'PDF', 'SHARP' and 'RESHARP' 
%       'refine'        : refine the RDF by 5th order polynomial fitting (all)
%       'tol'           : error tolerance (LBV or PDF)
%       'depth'         : multigrid level (LBV)
%       'peel'          : thickness of the boundary layer to be peeled off
%                         (LBV)
%       'b0dir'         : B0 direction (PDF)
%       'iteration'     : no. of maximum iteration for cgsolver (PDF)
%       'CGsolver'      : true for default cgsolver; false for matlab pcg
%                         solver (PDF)
%       'noisestd'      : noise standard deviation on the field map (PDF)                         
%       'radius'        : radius of the spherical mean value operation
%                         (SHARP or RESHARP)
%       'threshold'     : threshold used in Truncated SVD (SHARP)
%       'alpha'         : regularizaiton parameter used in Tikhonov
%                         (RESHARP)
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
                    [method, B0_dir, tol, iteration, CGdefault, N_std, refine] = parse_vararginPDF(varargin);
                    break
                case 'sharp'
                    [method, radius, threshold, refine] = parse_vararginSHARP(varargin);
                    break
                case 'resharp'
                    [method, radius, alpha, refine] = parse_vararginRESHARP(varargin);
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
        RDF = LBV(bkgField, mask, matrixSize, voxelSize, tol,depth,peel);
    case 'PDF'
        RDF = PDF(bkgField, N_std, mask,matrixSize,voxelSize, B0_dir, ...
            'tol', tol, 'iteration', iteration, 'CGsolver', CGdefault);
    case 'SHARP'
        RDF = SHARP(bkgField, mask, matrixSize, voxelSize, radius,threshold);
    case 'RESHARP'
        RDF = RESHARP(bkgField, mask, matrixSize, voxelSize, radius, alpha);
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
function [method, B0_dir, tol, iteration, CGdefault, N_std, refine] = parse_vararginPDF(arg)
method = 'PDF';
B0_dir = [0,0,1];
tol = 0.1;
iteration = 30;
CGdefault = true;
N_std = ones(size(bkgField))*1e-4;
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
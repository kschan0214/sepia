%% RDF = BackgroundRemovalMacro(bkgField,mask,matrixSize,voxelSize,varargin)
%
% Usage: 
%       RDF = BackgroundRemovalMacro(fieldMap,mask,matrixSize,voxelSize,...
%               'method','LBV','refine',true,'tol',1e-4,'depth',4,'peel',2);
%       RDF = BackgroundRemovalMacro(fieldMap,mask,matrixSize,voxelSize,...
%               'method','PDF','refine',true,'tol',0.1,'b0dir',[0,0,1],...
%               'iteration',100,'CGsolver',true,'noisestd',N_std);
%       RDF = BackgroundRemovalMacro(fieldMap,mask,matrixSize,voxelSize,...
%               'method','SHARP','refine',true,'radius',4,'threshold',0.03);
%       RDF = BackgroundRemovalMacro(fieldMap,mask,matrixSize,voxelSize,...
%               'method','RESHARP','refine',true,'radius',4,'alpha',0.01);
%       RDF = BackgroundRemovalMacro(fieldMap,mask,matrixSize,voxelSize,...
%               'method','VSHARP','refine',true); 
%       RDF = BackgroundRemovalMacro(fieldMap,mask,matrixSize,voxelSize,...
%               'method','iHARPERELLA','refine',true,'iteration',100);
%
% Description: Wrapper for background field removal (default using LBV)
%   Flags:
%       'method'        : background revomal method, 
%                          'LBV', 'PDF', 'SHARP', 'RESHARP', 'VSHARPSTI, and
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
%       VSHARPSTI 
%       ----------------
%
%       VSHARPSTI 
%       ----------------
%       'radius'        : radius of sphere
%
%       iHARPERELLA
%       ----------------
%       'iteration'     : no. of maximum iterations
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 28 June 2017
% Date last modified: 8 September 2017
%
function RDF = BackgroundRemovalMacro(totalField,mask,matrixSize,voxelSize,varargin)
%% Parsing argument input flags
if ~isempty(varargin)
    for kvar = 1:length(varargin)
        if strcmpi(varargin{kvar},'method')
            switch lower(varargin{kvar+1})
                case 'lbv'
                    method = 'LBV';
                    [tol, depth, peel, refine] = parse_varargin_LBV(varargin);
                    break
                case 'pdf'
                    method = 'PDF';
                    [B0_dir, tol, iteration, CGdefault, N_std, refine] = parse_varargin_PDF(varargin);
                    if isempty(N_std)
                        N_std = ones(matrixSize)*1e-4;
                    end
                    break
                case 'sharp'
                    method = 'SHARP';
                    [radius, threshold, refine] = parse_varargin_SHARP(varargin);
                    break
                case 'resharp'
                    method = 'RESHARP';
                    [radius, alpha, refine] = parse_varargin_RESHARP(varargin);
                    break
                case 'vsharpsti'
                    method = 'VSHARPSTI';
                    [refine] = parse_varargin_VSHARPSTI(varargin);
                    break
                case 'iharperella'
                    method = 'iHARPERELLA';
                    [iteration, refine] = parse_varargin_iHARPERELLA(varargin);
                    break
                case 'vsharp'
                    method = 'VSHARP';
                    [radius, refine] = parse_varargin_VSHARP(varargin);
            end
        end
    end
else
    % predefine paramater: if no varargin, use LBV
    disp('No method selected. Using the default setting:');
    method = 'LBV';
    tol = 1e-4;
    depth = 4;
    peel = 1;
    refine = false;
end

disp(['The following method is being used: ' method]);

%% background field removal
switch method
    case 'LBV'
        RDF = LBV(totalField,mask,matrixSize,voxelSize,tol,depth,peel);
    case 'PDF'
        RDF = PDF(totalField,mask,matrixSize,voxelSize,'b0dir',B0_dir,...
            'tol', tol,'iteration', iteration,'CGsolver', CGdefault,'noisestd',N_std);
    case 'SHARP'
        RDF = SHARP(totalField, mask, matrixSize, voxelSize, radius,threshold);
    case 'RESHARP'
        RDF = RESHARP(totalField, mask, matrixSize, voxelSize, radius, alpha);
    case 'VSHARPSTI'
        RDF = V_SHARP(totalField, mask,'voxelsize',double(voxelSize(:))');
    case 'iHARPERELLA'
        RDF = iHARPERELLA(totalField, mask,'voxelsize',voxelSize,'niter',iteration);
    case 'VSHARP'
        [RDF,~] = BKGRemovalVSHARP(totalField,mask,voxelSize,'radius',radius);
end

%% If refine is needed, do it now
if refine
     [~,RDF,~]=PolyFit(RDF,RDF~=0,5);
end

end

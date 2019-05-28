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
%       'refine'        : refine the RDF by 4th order polynomial fitting
%       'erode'         : number of voxels to be erode from edges
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
% Date modified: 29 September 2017
% Date modified: 1 April 2019
% Date modified: 24 May 2019
%
function RDF = BackgroundRemovalMacro(totalField,mask,matrixSize,voxelSize,varargin)
matrixSize = matrixSize(:).';
voxelSize = voxelSize(:).';

%% default algorithm parameters
refine          = false;
erode_radius    = 0;

%% Parsing argument input flags
if ~isempty(varargin)
    for kvar = 1:length(varargin)
        if strcmpi(varargin{kvar},'method')
            switch lower(varargin{kvar+1})
                case 'lbv'
                    method = 'LBV';
                    [tol, depth, peel] = parse_varargin_LBV(varargin);
%                     break
                case 'pdf'
                    method = 'PDF';
%                     [B0_dir, tol, iteration, CGdefault, N_std, refine] = parse_varargin_PDF(varargin);
                    [B0_dir, tol, iteration, padSize, N_std] = parse_varargin_PDF(varargin);
                    if isempty(N_std)
                        N_std = ones(matrixSize)*1e-4;
                    end
%                     break
                case 'sharp'
                    method = 'SHARP';
                    [radius, threshold] = parse_varargin_SHARP(varargin);
%                     break
                case 'resharp'
                    method = 'RESHARP';
                    [radius, alpha] = parse_varargin_RESHARP(varargin);
%                     break
                case 'vsharpsti'
                    method = 'VSHARPSTISuite';
                    [radius] = parse_varargin_VSHARPSTI(varargin);
%                     break
                case 'iharperella'
                    method = 'iHARPERELLA';
                    [iteration] = parse_varargin_iHARPERELLA(varargin);
%                     break
                case 'vsharp'
                    method = 'VSHARP';
                    [radius] = parse_varargin_VSHARP(varargin);
            end
        end
        if strcmpi(varargin{kvar},'refine')
            refine = varargin{kvar+1};
        end
        if strcmpi(varargin{kvar},'erode')
            erode_radius = varargin{kvar+1};
        end
    end
else
    % predefine paramater: if no varargin, use LBV
    disp('No method selected. Using the default setting:');
    method = 'LBV';
    tol = 1e-4;
    depth = 4;
    peel = 1;
end

% use single precision to reduce memory usage
totalField  = single(totalField);
mask       	= single(mask);
voxelSize   = single(voxelSize);
TE          = single(TE);
matrixSize  = single(matrixSize);
if exist('N_std','var')
    N_std = single(N_std);
end

disp(['The following method is being used: ' method]);

%% background field removal
% add path
sepia_addpath(method);

% zero padding for odd number dimension
totalField  = zeropad_odd_dimension(totalField,'pre');
mask        = zeropad_odd_dimension(mask,'pre');
if exist('N_std','var')
    N_std   = zeropad_odd_dimension(N_std,'pre');
end
matrixSize_new = size(totalField);
matrixSize_new = single(matrixSize_new);


disp('The following parameters are being used...');

switch method
    case 'LBV'
        disp(['Tolerance = ' num2str(tol)]);
        disp(['Depth = ' num2str(depth)]);
        disp(['Peel = ' num2str(peel)]);
        RDF = LBV(totalField,mask,matrixSize_new,voxelSize,tol,depth,peel);
        deleteme = dir('mask*.bin');
        delete(deleteme(1).name);
%         system(['rm ' deleteme.folder filesep deleteme.name]);
    case 'PDF'
        disp(['Tolerance = ' num2str(tol)]);
        disp(['Maximum iterations = ' num2str(iteration)]);
%         disp(['CGsolver = ' num2str(CGdefault)]);
%         RDF = PDF(totalField,mask,matrixSize,voxelSize,'b0dir',B0_dir,...
%             'tol', tol,'iteration', iteration,'CGsolver', CGdefault,'noisestd',N_std);
        RDF = PDF(totalField,N_std,mask,matrixSize_new,voxelSize,B0_dir,tol,...
            iteration,'imagespace',padSize);
    case 'SHARP'
        disp(['Radius(voxel) = ' num2str(radius)]);
        disp(['Threshold = ' num2str(threshold)]);
        RDF = SHARP(totalField, mask, matrixSize_new, voxelSize, radius,threshold);
    case 'RESHARP'
        disp(['Radius(voxel) = ' num2str(radius)]);
        disp(['Lambda = ' num2str(alpha)]);
        RDF = RESHARP(totalField, mask, matrixSize_new, voxelSize, radius, alpha);
        mask_RDF = SMV(mask, matrixSize_new, voxelSize, radius)>0.999;
        RDF = RDF .* mask_RDF;
    case 'VSHARPSTISuite'
        disp(['SMV size (mm): ' num2str(radius)]);
        RDF = V_SHARP(totalField, mask,'voxelsize',single(voxelSize(:))','smvsize',radius);
    case 'iHARPERELLA'
        disp(['Maximum iterations = ' num2str(iteration)]);
        RDF = iHARPERELLA(totalField, mask,'voxelsize',voxelSize,'niter',iteration);
    case 'VSHARP'
        disp(['Radius range(voxel) = ' num2str(radius)]);
        [RDF,~] = BKGRemovalVSHARP(totalField,mask,matrixSize_new,'radius',radius);
end

%% If refine is needed, do it now
if refine
    fprintf('Performing polynomial fitting...');
    % PolyFit required data to be double type
    RDF = double(RDF);
    [~,RDF,~]=PolyFit(RDF,RDF~=0,4);
    RDF = single(RDF);
    fprintf('Done!\n')
end

% get non-zero mask
if erode_radius > 0
    fprintf(['Eroding ' num2str(erode_radius) ' voxel(s) from edges...']);
    maskFinal = RDF ~=0;
    maskFinal = imfill(maskFinal,'holes');
    maskFinal = imerode(maskFinal,strel('sphere',erode_radius));
    % also remove the mask on the edges
    maskFinal(:,:,end-erode_radius:end) = 0;
    maskFinal(:,:,1:erode_radius)       = 0;
    maskFinal(:,end-erode_radius:end,:) = 0;
    maskFinal(:,1:erode_radius,:)       = 0;
    maskFinal(end-erode_radius:end,:,:) = 0;
    maskFinal(1:erode_radius,:,:)       = 0;
    RDF = RDF .* single(maskFinal);
    fprintf('Done!\n')
end

% remove zero padding 
RDF = zeropad_odd_dimension(RDF,'post',matrixSize);
% ensure the output is single to reduce memory usage
RDF = single(RDF);


end

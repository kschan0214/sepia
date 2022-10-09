%% [chi] = Wrapper_QSM_FANSI(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
%
% Input
% --------------
% localField    : local field map (tissue fields), in Hz
% mask          : signal mask
% matrixSize    : size of the input image
% voxelSize     : spatial resolution of each dimension of the data, in mm
% algorParam    : structure contains fields with algorithm-specific parameter(s)
% headerAndExtraData : structure contains extra header info/data for the algorithm
%
% Output
% --------------
% chi           : magnetic susceptibility map, in ppm
%
% Description: This is a wrapper function to access FANSI for SEPIA
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 8 March 2020
% Date modified: 20 Feb 2022 (v1.0)
%
%
function [chi] = Wrapper_QSM_FANSI(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
sepia_universal_variables;

% get algorithm parameters
algorParam = check_and_set_algorithm_default(algorParam);
method     = algorParam.qsm.method;
options.tol_update      = algorParam.qsm.tol;
options.maxOuterIter    = algorParam.qsm.maxiter;
options.mu2             = algorParam.qsm.mu2;
options.isWeakHarmonic  = algorParam.qsm.isWeakHarmonic;
options.beta            = algorParam.qsm.beta;
options.muh             = algorParam.qsm.muh;
alpha1                  = algorParam.qsm.lambda;
mu1                     = algorParam.qsm.mu1;

% need further decision
gradient_mode   = algorParam.qsm.gradient_mode;
constraint     	= algorParam.qsm.constraint;
solver          = algorParam.qsm.solver;

% get extra data such as magnitude/weights/B0 direction/TE/etc.
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
b0dir = headerAndExtraData.sepia_header.B0_dir;
b0    = headerAndExtraData.sepia_header.B0;
wmap  = get_variable_from_headerAndExtraData(headerAndExtraData, 'weights', matrixSize);
magn  = get_variable_from_headerAndExtraData(headerAndExtraData, 'magnitude', matrixSize);

% add path
sepia_addpath('FANSI');
addpath(fullfile(SEPIA_HOME,'misc','qsm_algorithm','FANSI'));

%% preparation
% algorithm parameters
switch lower(gradient_mode)
    case 'vector field'
        options.gradient_mode = 0;
    case 'l1 norm'
        options.gradient_mode = 1;
    case 'l2 norm'
        options.gradient_mode = 2;
    case 'none'
        % no gradient_mode
end

if strcmpi(constraint, 'tv')
    options.tgv = false;    % tv
else
    options.tgv = true;     % tgv
end

if strcmpi(solver,'linear')
    options.nonlinear = false;  % linear
else
    options.nonlinear = true;   % non-linear
end

% data 
% if both data are loaded
if ~isempty(magn) && ~isempty(wmap)
    disp('Both weighting map and magnitude images are loaded.');
    disp('Only the weighing map will be used.');
    
end
% if only magnitude images are loaded
if ~isempty(magn) && isempty(wmap)
    disp('The normalised RMS in time dimension of magnitude image will be used as the weighting map.');
    
    tmp     = sqrt(mean(magn.^2,4));
    wmap    = (tmp./max(tmp(:))) .* (mask); 
    
    clear tmp
end
% if nothing is loaded
if isempty(magn) && isempty(wmap)
    warning('Providing a weighing map or magnitude images can potentially improve the QSM map quality.');
    wmap = mask;
    options = rmfield(options, 'gradient_mode');
end

% masking weights
% wmap_min = min(wmap(wmap>0));
% wmap(wmap<wmap_min) = wmap_min;
wmap = wmap.*mask;

if ~isempty(magn)
    clear magn
end

%% Display algorithm parameters
disp('The following parameters are being used...');
disp(['Tolerance                    = ' num2str(options.tol_update)]);      
disp(['Max. iteration               = ' num2str(options.maxOuterIter)]);    
disp(['Fidelity consistancy (mu2)	= ' num2str(options.mu2)]);
disp(['Gradient L1 penalty (alpha1) = ' num2str(alpha1)]);
disp(['Gradient consistancy	(mu)    = ' num2str(mu1)]);          
% if isfield(options,'gradient_mode'); disp(['Gradient mode           = ' options.gradient_mode]); end
disp(['Gradient mode                = ' gradient_mode]);
disp(['Constraint                   = ' algorParam.qsm.constraint]);
disp(['Solver                       = ' algorParam.qsm.solver]);

disp(['Use weak-harmonic field regularisation?	= ' num2str(options.isWeakHarmonic)]);
if options.isWeakHarmonic
    disp(['Harmonic constraint          = ' num2str(options.beta)]);
    disp(['Harmonic consistancy         = ' num2str(options.muh)]);
end

%% main
% FANSI default parameters are optimised for ppm
localField = localField/(b0*gyro);

SNR     = 1000;     % Assume high SNR, mainly for avoiding division by zeros if regweight is used
noise   = 1/SNR;    % mainly to avoid division of zero
try
    % FANSI v1
%     noise   = 0;
    chi = FANSI_4sepia_v1(localField,wmap,voxelSize,alpha1,mu1,noise,options,b0dir);
    
catch
    % TODO 20220217: v2 need further testing!!
    
    % FANSI v3: scale of local field makes substantial differences with TGV
    % constraint
    % convert ppm to radian, similar to the demo scripts of FANSI
    phase_scale = 2*pi*gyro*b0 * headerAndExtraData.sepia_header.TE(end);
    localField = localField * phase_scale;
    
    % account for variable name changes
    if isfield(options,'gradient_mode'); options.gradientMode    = options.gradient_mode; end
    options.update          = options.tol_update;       options = rmfield(options, 'tol_update');       
    options.iterations      = options.maxOuterIter;     options = rmfield(options, 'maxOuterIter'); 
    options.voxelSize       = voxelSize;
    options.B0_dir          = b0dir;
    options.mu            	= mu1;       % updated v2
    options.noise           = noise; 
    options.isNonlinear     = options.nonlinear; 
    options.isTGV           = options.tgv; 
    % options.kernelMode   % SEPIA does not provide option of kernalMode    
    
    chi = FANSI_4sepia_v3( localField, wmap, alpha1, options );
    % chi = FANSI_4sepia(localField,wmap.*mask,voxelSize,alpha1,noise,options,b0dir);
    chi = chi.x;
    chi = chi / phase_scale; % radian to ppm
    
end

% chi = FANSI( localField, wmap, voxelSize, alpha1, noise, options, b0dir );

chi = chi .* mask;


end

%% set default parameter if not specified
function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;

try algorParam2.qsm.tol             = algorParam.qsm.tol;           catch; algorParam2.qsm.tol              = 0.1; end
try algorParam2.qsm.maxiter         = algorParam.qsm.maxiter;       catch; algorParam2.qsm.maxiter          = 150; end
try algorParam2.qsm.mu2             = algorParam.qsm.mu2;           catch; algorParam2.qsm.mu2              = 1; end
try algorParam2.qsm.lambda          = algorParam.qsm.lambda;        catch; algorParam2.qsm.lambda           = 4e-4; end
try algorParam2.qsm.mu1             = algorParam.qsm.mu1;           catch; algorParam2.qsm.mu1              = 4e-2; end
try algorParam2.qsm.solver          = algorParam.qsm.solver;        catch; algorParam2.qsm.solver           = 'Non-linear'; end
try algorParam2.qsm.constraint      = algorParam.qsm.constraint;	catch; algorParam2.qsm.constraint       = 'tv'; end
try algorParam2.qsm.gradient_mode   = algorParam.qsm.gradient_mode;	catch; algorParam2.qsm.gradient_mode	= 'Vector field'; end
try algorParam2.qsm.isWeakHarmonic  = algorParam.qsm.isWeakHarmonic;catch; algorParam2.qsm.isWeakHarmonic	= 0; end
try algorParam2.qsm.beta            = algorParam.qsm.beta;          catch; algorParam2.qsm.beta             = 150; end
try algorParam2.qsm.muh             = algorParam.qsm.muh;           catch; algorParam2.qsm.muh              = 3; end

end
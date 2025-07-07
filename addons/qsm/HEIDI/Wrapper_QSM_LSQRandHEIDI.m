%% [chi] = Wrapper_QSM_LSQRandHEIDI(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
%
% Input
% --------------
% localField    : local field map (tissue fields), in Hz
% mask          : signal mask
% matrixSize    : size of the input image
% voxelSize     : spatial resolution of each dimension of the data, in mm
% algorParam    : structure contains fields with algorithm-specific parameter(s)
% headerAndExtraData : structure contains extra header info/data for the algorithm
%                      MUST contain magnitude image, B0 direction (1x3),
%                      magnetic field strength (B0), echo times(ms; 1xn vector), 
%                      and rotation matrix (3x3).
%
% Output
% --------------
% chi           : magnetic susceptibility map, in ppm
%
% Description: This is a wrapper function to access HEIDI for SEPIA
%
% Created by: Fahad Salman & Ferdinand Schweser @ University at Buffalo
% Edited by: Kwok-Shing Chan @ MGH
% kchan2@mgh.harvard.edu
% Date created: 14 June 2025
% Date modified: 6 July 2025
%
%
function [chi] = Wrapper_QSM_LSQRandHEIDI(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)  % 20250614 KC: just make sure the input arguments have the same default input by design

sepia_universal_variables;

% 20250614 KC: system check
if ~isunix || ismac; error('HEIDI add-on on SEPIA currently supports Linux systems only'); end  % check if it is linux system
HEIDI_HOME = fullfile(SEPIA_HOME,'..','external','HEIDI_SEPIAready');                                      % temporary hard-coded, will update

% 20250614 KC: check whether we have execute right to GradientAnisotropicDiffusionImageFilter if not then give the right
shellcommand = [ 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:' fullfile(HEIDI_HOME,'HEIDI') ';' ...
fullfile(HEIDI_HOME,'HEIDI','GradientAnisotropicDiffusionImageFilter')];
status = system(shellcommand);
if status == 126    % if no executable right then grant the right
    shellcommand = [ 'chmod a+x ' fullfile(HEIDI_HOME,'HEIDI','GradientAnisotropicDiffusionImageFilter')];
    status = system(shellcommand);
end

%% Preprocessing for HEIDI using LSQR begins here

% get extra data such as magnitude/weights/B0 direction/TE/etc.
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);

magn = get_variable_from_headerAndExtraData(headerAndExtraData, 'magnitude', matrixSize);
Mask_CSF = [];

if ~isempty(magn) % maybe add some statement regarding multi-echo data otherwise this wont be possible

    if size(magn,4) == 1
        error('For magnitude data, only multi-echo data is supported. Please select a 4D [x,y,z,t] magnitude dataset or remove it from the input name');
    end
    disp('Extracting CSF mask....');
    % R2* mapping
    % r2s         = arlo(TE,iMag);
    r2s         = arlo(headerAndExtraData.sepia_header.TE,magn);    % 20250628: correct names
    Mask_CSF    = extract_CSF(r2s,mask,voxelSize)>0;

    clear r2s
end

% obtain header info
options.magneticFieldStrength   = headerAndExtraData.sepia_header.B0;
options.echoTime                = headerAndExtraData.sepia_header.TE *1e3; % 20250614 KC: HEIDI expect TE in ms
options.voxelAspectRatio        = headerAndExtraData.sepia_header.voxelSize;
options.B0_dir                  = headerAndExtraData.sepia_header.B0_dir;

% check and set default algorithm parameters
algorParam                                      = check_and_set_algorithm_default(algorParam,options);

% get algorithm parameters for LSQR
options.tolerance                               = algorParam.qsm.tolerance;
options.maxiter                                 = algorParam.qsm.maxiter;       % 20250614 KC: I suggest change this to 'maxiter' to make it consistent across all QSM methods available on SEPIA
options.offsetUseBool                           = algorParam.qsm.offsetUseBool;
options.isFourierDomainFormula                  = algorParam.qsm.isFourierDomainFormula;
options.TikhonovRegularizationSusceptibility    = algorParam.qsm.TikhonovRegularizationSusceptibility;
options.solvingType                             = algorParam.qsm.solvingType;
options.DipoleFilter                            = algorParam.qsm.DipoleFilter;
options.residualWeighting                       = algorParam.qsm.residualWeighting;

% 20250614 KC: don't need to add path for addon
% add path
addpath(fullfile(HEIDI_HOME,'LSQR'));

% 20250614 KC: magn need to be 3D, compress multu-echo data as MEDI
magn = sqrt(sum(abs(magn).^2,4));

initializeswitoolbox('docker',true)

mids;

DEFAULT_ANTIALIASINGFRAMEWIDTH = 40;
% mask
[mask,zcInfo] = zealouscrop(mask);
[mask,pcInfo] = prepareforconvolution(mask,[],DEFAULT_ANTIALIASINGFRAMEWIDTH);
% weight
spatialWeightingMatrix = zealouscrop(magn,zcInfo);
spatialWeightingMatrix = prepareforconvolution(spatialWeightingMatrix,[],DEFAULT_ANTIALIASINGFRAMEWIDTH);
% local feild
localField = zealouscrop(localField,zcInfo);
localField = prepareforconvolution(localField,[],DEFAULT_ANTIALIASINGFRAMEWIDTH);
% CSF mask
if ~isempty(Mask_CSF)
    Constraint.residual = zealouscrop(Mask_CSF,zcInfo);
    Constraint.residual = prepareforconvolution(Constraint.residual,[],DEFAULT_ANTIALIASINGFRAMEWIDTH);
else
    Constraint = [];
end

% 20250614: convert local field unit for HEIDI (expects in phase)
localField = localField /(2*pi);

% main - LSQR
[chi_raw,~, ~,~,~,~,~,~,~,phaseDiscrepancy] = computesusceptibility_lsqr(localField,options,logical(mask),Constraint,spatialWeightingMatrix);

%% Now HEIDI begins
% 20250614: don't need to add path for addon
rmpath(fullfile(HEIDI_HOME,'LSQR'));
addpath(fullfile(HEIDI_HOME,'HEIDI'));

initializeswitoolbox('docker',true)

mids;

% Although most of the options are similar, just to have a clear direction
% in terms of the options being used specifically for HEIDI, clear all 
% options to ensure the options are HEIDI-specific
clear options

% obtain header info
options.magneticFieldStrength   = headerAndExtraData.sepia_header.B0;
options.echoTime                = headerAndExtraData.sepia_header.TE(end) *1e3; % 20250614 KC: HEIDI expect TE in ms
options.voxelAspectRatio        = headerAndExtraData.sepia_header.voxelSize;
options.B0_dir                  = headerAndExtraData.sepia_header.B0_dir;

% get algorithm parameters
options.isFourierDomainFormula                  = algorParam.qsm.isFourierDomainFormula;
options.TikhonovRegularizationSusceptibility    = algorParam.qsm.TikhonovRegularizationSusceptibility;
options.solvingType                             = algorParam.qsm.solvingType;
options.DipoleFilter                            = algorParam.qsm.DipoleFilter;

options.useHEIDI                                = true;
options.PostProcCone.threshold                  = algorParam.qsm.PostProcCone_threshold;
options.PostProcCone.tol                        = algorParam.qsm.PostProcCone_tol;
options.PostProcCone.tolEnergy                  = algorParam.qsm.PostProcCone_tolEnergy; % makes much better results

% main - HEIDI
[chi,~, ~,~,~,~,~,~,~,phaseDiscrepancy] = computesusceptibility_heidi(chi_raw,localField,options,logical(mask),Constraint,spatialWeightingMatrix);

phase2susceptibilityFactor=phase2susceptibilityfactor([], options);
chi = (chi .* phase2susceptibilityFactor * 1e6) .* mask; 

chi = prepareforconvolution(chi,[],[],pcInfo);
chi = zealouscrop(chi,zcInfo,[],true);

end

%% set default parameter if not specified
function algorParam2 = check_and_set_algorithm_default(algorParam,options)

    algorParam2 = algorParam;
    
    try algorParam2.qsm.tolerance = algorParam.qsm.tolerance; % LSQR
        catch; algorParam2.qsm.tolerance = 1e-5; end
    try algorParam2.qsm.maxiter = algorParam.qsm.maxiter; % LSQR
        catch; algorParam2.qsm.maxiter = 400; end
    try algorParam2.qsm.offsetUseBool = algorParam.qsm.offsetUseBool; % LSQR
        catch; algorParam2.qsm.offsetUseBool = true; end
    try algorParam2.qsm.isFourierDomainFormula = algorParam.qsm.isFourierDomainFormula; % LSQR+HEIDI
        catch; algorParam2.qsm.isFourierDomainFormula = false; end
    try algorParam2.qsm.TikhonovRegularizationSusceptibility = algorParam.qsm.TikhonovRegularizationSusceptibility; % LSQR+HEIDI
        catch; algorParam2.qsm.TikhonovRegularizationSusceptibility = []; end
        if strcmpi(algorParam2.qsm.TikhonovRegularizationSusceptibility,'default'); algorParam2.qsm.TikhonovRegularizationSusceptibility = []; end  % 20250614 KC: compatible to GUI
    try algorParam2.qsm.solvingType = algorParam.qsm.solvingType; % LSQR+HEIDI
        if strcmpi(algorParam2.qsm.solvingType,'default'); algorParam2.qsm.solvingType = []; end % 20250614 KC: compatible to GUI
        catch; algorParam2.qsm.solvingType = []; end
    try algorParam2.qsm.DipoleFilter = algorParam.qsm.DipoleFilter; % LSQR+HEIDI
        if strcmpi(algorParam2.qsm.DipoleFilter,'default'); algorParam2.qsm.DipoleFilter = []; end % 20250614 KC: compatible to GUI
        catch; algorParam2.qsm.DipoleFilter = []; end
    try algorParam2.qsm.residualWeighting = algorParam.qsm.residualWeighting; % LSQR
        catch; algorParam2.qsm.residualWeighting = 0.2 / 9.4 * options.magneticFieldStrength; end

    % 20250706 KC: for loading the config file back, to the GUI, only one level of structure under qsm is supported, so I updated the following variables 
    try algorParam2.qsm.PostProcCone_threshold = algorParam.qsm.PostProcCone_threshold; % HEIDI
        catch; algorParam2.qsm.PostProcCone_threshold = 0.1; end
    try algorParam2.qsm.PostProcCone_tol = algorParam2.qsm.PostProcCone_tol; % HEIDI
        catch; algorParam2.qsm.PostProcCone_tol = eps; end
    try algorParam2.qsm.PostProcCone_tolEnergy = algorParam2.qsm.PostProcCone_tolEnergy; % HEIDI
        catch; algorParam2.qsm.PostProcCone_tolEnergy = eps; end

end

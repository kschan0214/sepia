%% Thanks for your interest in our STI SUITE software package!
% This file contains templates to use the STI SUITE functions.
% Author: Wei Li, PhD, and Chunlei, PhD


%% [1] The Laplacian-based phase unwrapping
% Inputs:
% rawphase: raw phase
% voxelsize: Spatial resolution
% padsize: size for padarray to increase numerical accuracy

[Unwrapped_Phase, Laplacian]=MRPhaseUnwrap(rawphase,'voxelsize',voxelsize,'padsize',padsize);
% If optional inputs are not assigned, the following default values will be used.
voxelsize=[1 1 1];
padsize=[12 12 12];


%% [2] HARPERELLA: Integrated phase unwrapping and background phase removal
% Inputs:
% Laplacian: tissue Laplacian
% BrainMask: Binary Brain Mask
% voxelsize: Spatial resolution, e.g. SpatialRes=[1 1 1];
% niter: Number of Iterations, e.g. nIter= 40;

[TissuePhase]=iHARPERELLA(Laplacian,BrainMask,'voxelsize',voxelsize,'niter',niter);
% If optional inputs are not assigned, the following default values will be used.
voxelsize=[1 1 1];
niter=100;
NewMask=BrainMask;

%% [3] V-SHARP: background phase removal
% Inputs:
% Laplacian:	Raw phase Laplacian.
% voxelsize: voxel size 
% BrainMask: Binary Brain Mask

[TissuePhase,NewMask]=V_SHARP(Unwrapped_Phase,BrainMask,'voxelsize',voxelsize);

%% [4] iLSQR: Quantative Susceptibility Mapping
% Inputs:
% TissuePhase: tissue phase
% SpatialRes: Spatial resolution
% padsize: size for padarray to increase numerical accuracy
% H: the field direction, e.g. H=[0 0 1];
% ninter: number of iterations. e.g. niter=30.

% ------------ usage 1 -----------------------
[Susceptibility]= QSM_iLSQR(TissuePhase,BrainMask,'H',H,'voxelsize',voxelsize,'padsize',padsize,'niter',niter,'TE',TE,'B0',B0);
% If optional inputs are not assigned, the following default values will be used.
H=[0 0 1];
voxelsize=[1 1 1];
padsize=[0 0 0];
B0=3;
TE=40; 

% ------------ usage 2 -----------------------
params.H=H;
params.voxelsize = voxelsize;
params.padsize = padsize;
params.niter = max_niter;
params.cropsize = cropsize;
params.TE = TE;
params.B0 = B0;
params.tol_step1 = tol_step1;
params.tol_step2 = tol_step2;
params.Kthreshold = Kthreshold;

[Susceptibility]= QSM_iLSQR(TissuePhase,NewMask,'params',params);



%% [5] Fast QSM
% Inputs:
% TissuePhase: tissue phase
% voxelsize: Spatial resolution
% padsize: size for padarray to increase numerical accuracy
% H: the field direction, e.g. H=[0 0 1];

[susceptibility_fs]=FastQSM(TissuePhase,NewMask,'H',H,'voxelsize',voxelsize,'padsize',padsize,'TE',TE,'B0',B0);
% If optional inputs are not assigned, the following default values will be used.
voxelsize=[1 1 1];
H=[0 0 1];
B0=3;
TE=1;
padsize=[0 0 0];

%% [6] Susceptibility Tensor Imaging
% Inputs: 
% Phase4D: 4D dataset the 4th dimention is orientation. e.g. 
%          size(Phase4D)=[128,128,128,n];                   % n directions
% H_Matrix: the H matrix
%          e.g. H= [  0.0072   -0.8902   -0.4556            % direction 1
%                     0.1121   -0.2320   -0.9662            % direction 2
%                     ........
%                    -0.1384    0.7923    0.5942];          % direction n
% B0: B0, e.g. B0=3;
% TE_ms: TE in the unit of ms, e.g TE=30.
% parallel_Flag: parallel computing flag. recommend 'off'.

[SusceptibilityTensor]=STI_Parfor(Phase4D,H_Matrix,B0,TE_ms,parallel_Flag);



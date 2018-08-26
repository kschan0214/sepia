%% Thanks for your interest in our STI SUITE software package!
% This file contains templates to use the STI SUITE functions.
% Author: Hongjiang Wei, PhD, and Chunlei, PhD
%% [0] coil combination functions
% read 4D dataset [x y z coils] and combine multicoil, return mag, phase of each echo
% this multiple coils
% Input
% img        -  4D array, [Nx, Ny, Nz, Ncoil]
% SpatialRes -  resolution, e.g. [1 1 2]
% Output
% magniCombined        -  magnitude image
% phaseCombined        -  phase image

[magniCombined, phaseCombined] = Combine_Coils(img,SpatialRes);

%% 

[data,voxelsize,matrix_size,frequency,delta_TE,TEs,B0_vector,B0] = Read_DICOM_HW(FilePath);
mag = abs(data);
phase = angle(data);
% read dicom images from a single folder
%
% Input
% dicom file path
% Output
% data        -  complex 4D data
% voxelsize
% matrix_size
% frequency: central frequency
% delta_TE: TE spacing
% TEs: multi-TEs 
% B0_vector: H vector, e.g. [0 0 1]
% B0: field strength, e.g. 3


%% [1] The Laplacian-based phase unwrapping
% Inputs:
% rawphase: 3D raw phase
% voxelsize: Spatial resolution
% padsize: size for padarray to increase numerical accuracy

[Unwrapped_Phase, Laplacian]=MRPhaseUnwrap(rawphase,'voxelsize',voxelsize,'padsize',padsize);
% If optional inputs are not assigned, the following default values will be used.
voxelsize=[1 1 1];
padsize=[12 12 12];

%% [2] 2D V-SHARP: background phase removal for 2D QSM based on 2D GRE EPI scan
% Inputs:
% Unwrapped_Phase: 3D unwrapped phase using Laplacian unwrapping
% BrainMask: brain mask 
% voxelsize: spatial resoluiton 
% padsize: size for padarray to increase numerical accuracy
% smvsize: filtering size, default value = 12
% Output:
% 3D TissuePhase

TissuePhase = V_SHARP_2d(Unwrapped_Phase,BrainMask,'voxelsize',voxelsize,'padsize',padsize,'smvsize',smvsize);
smvsize=12;
voxelsize=[1 1 1];

%% [4] --------- iHARPERELLA using wrapped phase ---------------
% Inputs:
% rawphase: 3D raw phase 
% BrainMask: brain mask 
% voxelsize: spatial resoluiton 
% padsize: size for padarray to increase numerical accuracy
% niter: iteration number
%
% Output:
% 3D Filtered TissuePhase

[TissuePhase]=iHARPERELLA(rawphase,BrainMask,'voxelsize',voxelsize,'padsize',padsize,'niter',niter);
% If optional inputs are not assigned, the following default values will be used.
niter=40;
voxelsize=[1 1 1];
padsize=[12 12 12];

%% [5] fast STAR-QSM (~30 sec): Quantative Susceptibility Mapping
% Inputs:
% TissuePhase: tissue phase
% SpatialRes: Spatial resolution
% padsize: size for padarray to increase numerical accuracy
% H: the field direction, e.g. H=[0 0 1];
 Susceptibility = QSM_star(TissuePhase,NewMask,'TE',TE,'B0',B0,'H',H,'padsize',padsize,'voxelsize',voxelsize);
% If optional inputs are not assigned, the following default values will be used.
H=[0 0 1];
voxelsize=[1 1 1];
padsize=[12 12 12];
B0=3;
TE=40; 

%% [6] iLSQR: Quantative Susceptibility Mapping
% Inputs:
% TissuePhase: tissue phase
% SpatialRes: Spatial resolution
% padsize: size for padarray to increase numerical accuracy
% H: the field direction, e.g. H=[0 0 1];
[Susceptibility]= QSM_iLSQR(TissuePhase,NewMask,'TE',TE,'B0',B0,'H',H,'padsize',padsize,'voxelsize',voxelsize);
% If optional inputs are not assigned, the following default values will be used.
H=[0 0 1];
voxelsize=[1 1 1];
padsize=[12 12 12];
B0=3;
TE=40; 

%% [67] Susceptibility Tensor Imaging
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



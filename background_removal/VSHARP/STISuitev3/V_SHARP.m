%  3D background phase removal for 3D GRE data, 
%  Inputs:
%  Unwrapped_Phase: 3D unwrapped phase using Laplacian unwrapping
%  BrainMask: brain mask 
%  voxelsize: spatial resoluiton 
%  padsize: size for padarray to increase numerical accuracy
%  smvsize: filtering size, default value = 12
% 
%  Output:
%  TissuePhase: 3D Filtered TissuePhase
%  NewMask: eroded mask
%  For example: 
%  [TissuePhase,NewMask]=V_SHARP(Unwrapped_Phase,BrainMask,'voxelsize',voxelsize,'smvsize',smvsize);
%  Hongjiang Wei, PhD
%  Chunlei Liu, PhD
%  University of California, Berkeley
%

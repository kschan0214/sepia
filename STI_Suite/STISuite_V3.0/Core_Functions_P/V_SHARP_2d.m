%  2D background phase removal for 2D GRE-EPI data, 
%  Inputs:
%  Unwrapped_Phase: 3D unwrapped phase using Laplacian unwrapping
%  BrainMask: brain mask 
%  voxelsize: spatial resoluiton 
%  padsize: size for padarray to increase numerical accuracy
%  smvsize: filtering size, default value = 12
%  Output:
%  3D TissuePhase
%  For example: 
%  TissuePhase = V_SHARP_2d(Unwrapped_Phase,BrainMask,'voxelsize',voxelsize,'padsize',padsize,'smvsize',smvsize);
%  Hongjiang Wei, PhD
%  Chunlei Liu, PhD
%  University of California, Berkeley
%

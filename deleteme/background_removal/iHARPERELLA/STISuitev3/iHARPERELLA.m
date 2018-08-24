%  intergrated phase unwrapping and background phase removal for 3D GRE data, 
%  Inputs:
%  rawphase: 3D raw phase 
%  BrainMask: brain mask 
%  voxelsize: spatial resoluiton 
%  padsize: size for padarray to increase numerical accuracy
%  niter: iteration number
% 
%  Output:
%  3D Filtered TissuePhase
%  For example: 
%  [TissuePhase]=iHARPERELLA(rawphase,BrainMask,'voxelsize',voxelsize,'padsize',padsize,'niter',niter);
%  Wei Li, PhD
%  Chunlei Liu, PhD
%  University of California, Berkeley
%  Wei Li, March 3,2014.
%

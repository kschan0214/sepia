%  star QSM for QSM calculation, 
%  Inputs:
%  TissuePhase: tissue phase
%  NewMask: brain mask
%  SpatialRes: Spatial resolution
%  padsize: size for padarray to increase numerical accuracy
%  H: the field direction, e.g. H=[0 0 1];
%  Output:
%  Susceptibility: QSM images
%  For example: 
%  Susceptibility = QSM_star(TissuePhase,NewMask,'TE',TE,'B0',B0,'H',H,'padsize',padsize,'voxelsize',voxelsize);
%  Hongjiang Wei, PhD
%  Chunlei Liu, PhD
%  University of California, Berkeley
%  Updated on 02/09/2014
%

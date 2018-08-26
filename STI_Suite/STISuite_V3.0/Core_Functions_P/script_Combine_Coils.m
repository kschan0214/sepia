%  Combine_Coils
%  % read 4D dataset [x y z coils] and combine multicoil, return mag, phase of each echo
% % Syntax
%   [magniCombined, phaseCombined] = Combine_Coils(img,SpatialRes)
% % Description
%  this multiple coils
%  Input
%  img        -  4D array, [Nx, Ny, Nz, Ncoil]
%  SpatialRes -  resolution, e.g. [1 1 2]
%  Output
%  magniCombined        -  magnitude image
%  phaseCombined        -  phase image
%  Chunlei Liu, PhD
%  05/26/2014
% %
% % parameters
% % parameters
%

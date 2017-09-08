%% function chi = qsmTKD(localField,mask,matrixSize,voxelSize,thre_tkd)
%
% Description: compute QSM based on threshol-based k-space division (TKD)
% Ref        : Wharton et al. MRM 63:1292-1304(2010)
%
% Input
% -----
%   localField      : local field perturbatios
%   mask            : user-defined mask
%   matrixSize      : image matrix size
%   voxelSize       : spatial resolution of image 
%   varargin        : flags with
%       'threshold'     -   threshold for k-space inversion 
%
% Output
% ______
%   chi             : QSM
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 24 March 2017
% Date last modified: 6 September 2017
%
function chi = qsmTKD(localField,mask,matrixSize,voxelSize,varargin)
%% Parsing varargin
thre_tkd = parse_varargin_TKD(varargin);

% display message
fprintf('Threshold for k-space division is %f \n',thre_tkd);

%% Core
% dipole kernel
kernel = DipoleKernel(matrixSize,voxelSize);

% initiate inverse kernel with zeros
kernel_inv = zeros(matrixSize);
% get the inverse only when value > threshold
kernel_inv( abs(kernel) > thre_tkd ) = 1 ./ kernel(abs(kernel) > thre_tkd);
% direct dipole inversion method
chi = real( ifftn( fftn(localField) .* kernel_inv ) ) .* mask;

end

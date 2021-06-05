%% chi = myQSM(localField,mask,matrixSize,voxelSize,thres,b0dir)
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
% Dates modified: 6 September 2017
%
function chi = myQSM(localField,mask,matrixSize,voxelSize,thres,b0dir)
% if the last two input variables are not provided, set to a default value
if nargin < 6
    b0dir = [0,0,1];
end
if nargin < 5
    thres = 0.15;
end

% display message
fprintf('Threshold for k-space division is %f \n',thres);

%% Core
% dipole kernel
kernel = DipoleKernel(matrixSize,voxelSize,b0dir);

% initiate inverse kernel with zeros
kernel_inv = zeros(matrixSize, 'like', matrixSize);

% get the inverse only when value > threshold
kernel_inv( abs(kernel) > thres ) = 1 ./ kernel(abs(kernel) > thres);

% direct dipole inversion method
chi = real( ifftn( fftn(localField) .* kernel_inv ) ) .* mask;

end

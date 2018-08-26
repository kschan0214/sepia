% Correction of the iField - removal of the echo-dependent linear phase
%   gradient to avoid non-2pi wrap-like artifacts in iFreq_raw 
%
%   [iField_corrected] = iField_correction(iField,voxel_size)
% 
%   output
%   iField_corrected - iField with echo-dependent linear gradient
%   removed from the phase
%
%   input
%   iField - the complex MR image
%   voxel_size - the size of a voxel
%
%   Created by Alexey V. Dimov and Pascal Spincemaille in 2017.01.13
%   Updated by Alexey V. Dimov on 2017.02.02

function [iField_corrected] = iField_correction_upd(iField,voxel_size)
matrix_size = size(squeeze(iField(:,:,:,1)));
Mask = genMask(iField,voxel_size);
pha = angle(iField);
mag = sqrt(sum(abs(iField).^2,4));

%Laplacian in echo direction
echo_laplacian = zeros([matrix_size, size(iField,4)-2]);
for i = 1 : size(echo_laplacian,4)
    echo_laplacian(:,:,:,i) = pha(:,:,:,i) - 2*pha(:,:,:,i+1) + pha(:,:,:,i+2);
end
echo_laplacian = cat(4, pha(:,:,:,2) - pha(:,:,:,1), echo_laplacian);

%Quick and dirty unwrapping
for i = 1 : size(echo_laplacian,4)
    echo_laplacian(:,:,:,i) = unwrapPhase(mag,squeeze(echo_laplacian(:,:,:,i)),matrix_size);
end

%Estimate parameters of the gradient
[g1, g2, in] = gffun(echo_laplacian,Mask);
[slope1, slope2, intrcp]  = fitgrad(g1, g2, in);

%Subtract artificially added phase
phasor = phasorprep(slope1, slope2, intrcp, matrix_size);
iField_corrected = iField.*exp(-1i*phasor);
end
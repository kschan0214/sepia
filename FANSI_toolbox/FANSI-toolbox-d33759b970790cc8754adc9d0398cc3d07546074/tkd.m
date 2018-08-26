% QSM through the Truncation of the Dipole Kernel method
%
% Based on the code by Bilgic Berkin at http://martinos.org/~berkin/software.html
% Last modified by Carlos Milovic in 2017.12.27
%

function chi_tkd = tkd( phase_use, mask_use, kernel, kthre, N, trunc )
%
% input:
% phase_use - local field map
% mask_use - binary 3D image that defines the ROI.
% kernel - dipole kernel in the frequency space
% kthre - threshold in the frequency space to truncate the kernel
% N - array size
%
% output:
% chi_tkd - susceptibility map.
%

if nargin < 6
    trunc = true;
end

if trunc == true
kernel_inv = zeros(N);
else
kernel_inv = ones(N)/kthre; 
end
    
kernel_inv( abs(kernel) > kthre ) = 1 ./ kernel(abs(kernel) > kthre);

%tic
    chi_tkd = real( ifftn( kernel_inv.* fftn(phase_use) ) ) .* mask_use; 
%toc


end

% QSM with a Tikhonov regularization (L2 norm of solution) method
%
% Based on the code by Bilgic Berkin at http://martinos.org/~berkin/software.html
% Last modified by Carlos Milovic in 2017.03.30
%

function chi_L2 = chiL2( phase_use, mask_use, kernel, beta, N )
%
% input:
% phase_use - local field map
% mask_use - binary 3D image that defines the ROI.
% kernel - dipole kernel in the frequency space
% beta - regularization parameter
% N - array size
%
% output:
% chi_L2 - susceptibility map.
%

[kx, ky, kz] = ndgrid(0:N(1)-1, 0:N(2)-1, 0:N(3)-1);
Ex = 1 - exp(2i .* pi .* kx / N(1));
Ey = 1 - exp(2i .* pi .* ky / N(2));
Ez = 1 - exp(2i .* pi .* kz / N(3));

Ext = conj(Ex);
Eyt = conj(Ey);
Ezt = conj(Ez);

E2 = Ext .* Ex + Eyt .* Ey + Ezt .* Ez;
K2 = abs(kernel).^2;


tic
    chi_L2 = real( ifftn(conj(kernel) .* fftn(phase_use) ./ (K2 + beta * E2)) ) .* mask_use;
toc


end

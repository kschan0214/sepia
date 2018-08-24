% Calculates the L1 norm of the difference of the gradients between two images
% 
% Last modified by Carlos Milovic in 2017.03.30
%
function [ mad ] = compute_mad_GX( chi_recon, chi_true )


    N = size(chi_true);
    [k1, k2, k3] = ndgrid(0:N(1)-1, 0:N(2)-1, 0:N(3)-1);
    E1 = 1 - exp(2i .* pi .* k1 / N(1));
    E2 = 1 - exp(2i .* pi .* k2 / N(2));
    E3 = 1 - exp(2i .* pi .* k3 / N(3)); 
    
    Gt1 = real(ifftn(E1.*fftn(chi_true)));
    Gt2 = real(ifftn(E2.*fftn(chi_true)));
    Gt3 = real(ifftn(E3.*fftn(chi_true)));
    Gr1 = real(ifftn(E1.*fftn(chi_recon)));
    Gr2 = real(ifftn(E2.*fftn(chi_recon)));
    Gr3 = real(ifftn(E3.*fftn(chi_recon)));
    
    d1 = Gt1 - Gr1;
    d2 = Gt2 - Gr2;
    d3 = Gt3 - Gr3;
    GX =  abs(d1) + abs(d2) + abs(d3);



mad = 100 * sum( GX(:) ) / sum( abs(Gt1(:)) + abs(Gt2(:)) + abs(Gt3(:)) );


end


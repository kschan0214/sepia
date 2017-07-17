% Calculates the L2 norm of the difference of the gradients between two images
% 
% Last modified by Carlos Milovic in 2017.03.30
%
function [ rmse ] = compute_rmse_GX( chi_recon, chi_true )


    N = size(chi_true);
    [k1, k2, k3] = ndgrid(0:N(1)-1, 0:N(2)-1, 0:N(3)-1);
    E1 = 1 - exp(2i .* pi .* k1 / N(1));
    E2 = 1 - exp(2i .* pi .* k2 / N(2));
    E3 = 1 - exp(2i .* pi .* k3 / N(3)); 
    
    Gt1 = real(ifftn(E1.*fftn(chi_true)));
    Gt2 = real(ifftn(E2.*fftn(chi_true)));
    Gr1 = real(ifftn(E1.*fftn(chi_recon)));
    Gr2 = real(ifftn(E2.*fftn(chi_recon)));
    
    d1 = Gt1 - Gr1;
    d2 = Gt2 - Gr2;
    GX = sqrt( (d1).^2 + (d2).^2 );



rmse = 100 * norm( GX(:) ) / norm(sqrt((Gt1(:)).^2 + (Gt2(:)).^2 ));


end


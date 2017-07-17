% Calculates the gradient of a 3D data cube, and outputs the gradient 
% magnitude or vector field
%
% Last modified by Carlos Milovic in 2017.03.30
%
function [ output ] = gradient_calc( m, data_mode )
%
% input:
% m - input image. Typically the magnitude data for regularization weighting in QSM.
% data_mode - 0 to return the vector field. 1 for the L1 norm, and 2 for the L2 norm.
%
% output:
% output - gradient image.
%
    N = size(m);

    [k1, k2, k3] = ndgrid(0:N(1)-1, 0:N(2)-1, 0:N(3)-1);
    E1 = 1 - exp(2i .* pi .* k1 / N(1));
    E2 = 1 - exp(2i .* pi .* k2 / N(2));
    E3 = 1 - exp(2i .* pi .* k3 / N(3)); 
    
    Gm1 = real(ifftn(E1.*fftn(m)));
    Gm2 = real(ifftn(E2.*fftn(m)));
    Gm3 = real(ifftn(E3.*fftn(m)));
    
    if data_mode == 0
        output = abs(Gm1); 
        output(:,:,:,2) = abs(Gm2);
        output(:,:,:,3) = abs(Gm3);
    elseif data_mode == 1
        output = abs( Gm1 ) + abs( Gm2 ) + abs( Gm3 );
    elseif data_mode == 2
        output = sqrt( Gm1.^2 + Gm2.^2 +Gm3.^2 );
    end
    
    
end


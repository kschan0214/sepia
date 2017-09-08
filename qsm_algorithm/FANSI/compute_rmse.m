% Calculates the L2 norm of the difference between two images
% 
% Last modified by Carlos Milovic in 2017.03.30
%
function [ rmse ] = compute_rmse( chi_recon, chi_true )


rmse = 100 * norm( chi_recon(:) - chi_true(:) ) / norm(chi_true(:));


end


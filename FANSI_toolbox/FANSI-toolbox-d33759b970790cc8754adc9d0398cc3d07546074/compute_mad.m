% Calculates the L1 norm of the difference between two images
% 
% Last modified by Carlos Milovic in 2017.03.30
%
function [ mad ] = compute_mad( chi_recon, chi_true )


mad = 100 * sum( abs(chi_recon(:) - chi_true(:)) ) / sum(abs((chi_true(:))));


end


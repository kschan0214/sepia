% Calculates several quality indices between two images.
% See README.txt for further information.
% Last modified by Carlos Milovic in 2017.03.30
%
function [ metrics ] = compute_metrics( chi_recon, chi_true )

metrics.rmse = compute_rmse( chi_recon, chi_true );
metrics.hfen = compute_hfen( chi_recon, chi_true );
metrics.ssim = compute_ssim( chi_recon, chi_true );
metrics.cc = compute_cc( chi_recon, chi_true );
metrics.mi = compute_mi( chi_recon, chi_true );
metrics.gxe = compute_rmse_GX( chi_recon, chi_true );


end


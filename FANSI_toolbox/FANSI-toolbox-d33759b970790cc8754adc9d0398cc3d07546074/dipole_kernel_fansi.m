% Continuous dipole frequency kernel for arbitrary orientation
function [ kernel ] = dipole_kernel_fansi( N, spatial_res, B0_dir )
%
% input:
% N - array size
% spatial_res - 
% mode - 0 for the continuous kernel proposed by Salomir, et al. 2003.
%        1 for the discrete kernel proposed by Milovic, et al. 2017.
%        2 for the Integrated Green function proposed by Jenkinson, et al. 2004
% B0_dir - main field direction, e.g. 0 0 1
%
% output:
% kernel - dipole kernel in the frequency space
%
% Created by Carlos Milovic, 30.03.2017
% Modified by Julio Acosta-Cabronero, 26.05.2017
% Last modified by Carlos Milovic, 27.12.2017

[ky,kx,kz] = meshgrid(-N(2)/2:N(2)/2-1, -N(1)/2:N(1)/2-1, -N(3)/2:N(3)/2-1);

kx = (kx / max(abs(kx(:)))) / spatial_res(1);
ky = (ky / max(abs(ky(:)))) / spatial_res(2);
kz = (kz / max(abs(kz(:)))) / spatial_res(3);

k2 = kx.^2 + ky.^2 + kz.^2;
k2(k2==0) = eps;


% JAC
kernel = fftshift( 1/3 - (kx*B0_dir(1) + ky*B0_dir(2) + kz*B0_dir(3)).^2 ./ k2 );    
kernel(1,1,1) = 0.0;

    

end


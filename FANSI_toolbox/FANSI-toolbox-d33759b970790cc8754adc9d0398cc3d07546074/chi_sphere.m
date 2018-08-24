% Calculates the susceptibility distribution and magnetization by a sphere.
% This generates output maps that show the sampled values at each voxel location.
%
% Last modified by Carlos Milovic in 2017.12.27
%
function [chi, Bnuc, Xs] = chi_sphere(FOV,N,center,radius, xin, xout)
%
% input:
% FOV - field of view in x, y, and z directions
% N   - no of samples in kx, ky, kz
% center - center index of the k-space (cx,cy,cx)
% radius - radius of the sphere in voxels.
% xin - susceptibility value inside the sphere
% xout - susceptibility value outside the sphere
%
% output:
% chi - sampled susceptibility distribution
% Bnuc - analytic nuclear magnetic field (using the Lorent's sphere approximation)
% Xs - analytic fourier transform of spherical distribution
%


kx = 1:N(1);
ky = 1:N(2);
kz = 1:N(3);

kx = kx - center(1);
ky = ky - center(2);
kz = kz - center(3);

delta_kx = FOV(1)/N(1);
delta_ky = FOV(2)/N(2);
delta_kz = FOV(3)/N(3);


kx = kx * delta_kx;
ky = ky * delta_ky;
kz = kz * delta_kz;

kx = reshape(kx,[length(kx),1,1]);
ky = reshape(ky,[1,length(ky),1]);
kz = reshape(kz,[1,1,length(kz)]);

kx = repmat(kx,[1,N(2),N(3)]);
ky = repmat(ky,[N(1),1,N(3)]);
kz = repmat(kz,[N(1),N(2),1]);

k2 = kx.^2 + ky.^2 + kz.^2;

chi = zeros(N);

chi(k2 > radius*radius) = xout;

chi(k2 <= radius*radius) = xin;

% Bnuc calculated as in Salomir, 2003.
dX = xin-xout;
Bnuc = (1-chi*2/3);
Bmac = (radius^3)*dX*(2*kz.*kz-kx.*kx-ky.*ky)./((3+dX)*k2.^(5/2));
Bmac(k2 <= radius*radius) = 0;%2*dX/(3+dX);
Bnuc = Bnuc.*Bmac;


q = radius*sqrt( (kx.^2)/N(1)^2 + (ky.^2)/N(2)^2 + (kz.^2)/N(3)^2);
q(k2==0) = 1e-6;
p = 2*pi*q;
Xs = -p.*cos(p)+sin(p);
Xs = xin*Xs./(2*pi*pi*q.^3);
%Xs( q == 0) = xin*4*pi/(3*radius^3);

Xs = radius^3*Xs/sqrt(N(1)*N(2)*N(3));

end

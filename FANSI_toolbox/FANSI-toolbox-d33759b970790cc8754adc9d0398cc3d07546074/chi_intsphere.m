% Calculates the susceptibility distribution and magnetization by a sphere.
% This generates output maps that show the integrated values inside the 
% defined voxels. This provides a more realistic analytic ground truth, 
% based on the acquisition model on MRI devises.
%
% Last modified by Carlos Milovic in 2017.03.30
%
function [chi, Bnuc] = chi_intsphere(FOV,N,center,radius, xin, xout)
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
% chi - per voxel integrated susceptibility distribution
% Bnuc - analytic nuclear magnetic field (using the Lorent's sphere approximation)
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

M = zeros(N);
Mi = zeros(N);
Mo = zeros(N);
ri = max(radius-2,0.0);
ro = radius+15;
Mi(k2 >= ri^2) = 1;
Mo(k2<=ro^2) = 1;
M = Mi.*Mo;

chi = zeros(N);
chi(k2 > radius*radius) = xout;
chi(k2 <= radius*radius) = xin;

dX = xin-xout;
Bnuc = zeros(N);
Bnuc = (1-chi*2/3);
Bmac = (radius^3)*dX*(2*kz.*kz-kx.*kx-ky.*ky)./((3+dX)*k2.^(5/2));
Bmac(k2 <= radius*radius) = 0;%2*dX/(3+dX);
Bnuc = Bnuc.*Bmac;

overs = 15; % oversample factor. Near the boundaries of the sphere use a grid to approximate
            % the integrated susceptibility and field by averaging analytic values inside a 
            % voxel. Away from the boundaries, the central value is considered the mean value
for i=1:N(1)
    for j=1:N(2)
        for k=1:N(3)
            
            if (M(i,j,k) == 1)
            
            x(1) = kx(i,j,k)-0.5*delta_kx;
            x(2) = kx(i,j,k)+0.5*delta_kx;
            y(1) = ky(i,j,k)-0.5*delta_ky;
            y(2) = ky(i,j,k)+0.5*delta_ky;
            z(1) = kz(i,j,k)-0.5*delta_kz;
            z(2) = kz(i,j,k)+0.5*delta_kz;
            
            
            dx = 1:overs;
            dy = 1:overs;
            dz = 1:overs;
            dx = (x(2)-x(1))*dx/overs +x(1);
            dy = (y(2)-y(1))*dy/overs +y(1);
            dz = (z(2)-z(1))*dz/overs +z(1);
            
dx = reshape(dx,[length(dx),1,1]);
dy = reshape(dy,[1,length(dy),1]);
dz = reshape(dz,[1,1,length(dz)]);

dx = repmat(dx,[1,overs,overs]);
dy = repmat(dy,[overs,1,overs]);
dz = repmat(dz,[overs,overs,1]);

d2 = dx.^2 + dy.^2 + dz.^2;
            
            dchi = zeros([overs overs overs]);
            
dchi(d2 > radius*radius) = xout;

dchi(d2 <= radius*radius) = xin;
            
            chi(i,j,k) = mean(dchi(:));
            
            dBmac = (radius^3)*dX*(2*dz.*dz-dx.*dx-dy.*dy)./((3+dX)*d2.^(5/2));
            dBmac(d2 <= radius*radius) = 0;
            dBnuc = (1-dchi*2/3);
            Bnuc(i,j,k) = mean(dBmac(:).*dBnuc(:));
            end
        end
    end
end



end

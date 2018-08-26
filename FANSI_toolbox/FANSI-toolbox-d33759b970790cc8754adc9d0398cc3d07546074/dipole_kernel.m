% Dipole frequency kernel for main fields in the z axis
function [ kernel ] = dipole_kernel( N, spatial_res, mode )
%
% input:
% N - array size
% spatial_res - 
% mode - 0 for the continuous kernel proposed by Salomir, et al. 2003.
%        1 for the discrete kernel proposed by Milovic, et al. 2017.
%        2 for the Integrated Green function proposed by Jenkinson, et al. 2004
%
% output:
% kernel - dipole kernel in the frequency space
%
% Last modified by Carlos Milovic, 27.12.2017

[ky,kx,kz] = meshgrid(-N(2)/2:N(2)/2-1, -N(1)/2:N(1)/2-1, -N(3)/2:N(3)/2-1);


if mode == 0 % Continuous kernel

kx = (kx / max(abs(kx(:)))) / spatial_res(1);
ky = (ky / max(abs(ky(:)))) / spatial_res(2);
kz = (kz / max(abs(kz(:)))) / spatial_res(3);

k2 = kx.^2 + ky.^2 + kz.^2;

R_tot = eye(3);

kernel = fftshift( 1/3 - (kx * R_tot(3,1) + ky * R_tot(3,2) + kz * R_tot(3,3)).^2 ./ (k2 + eps) ); 
%kernel13 = fftshift( kz.* kx  ./ (k2 + eps) );      
%kernel23 = fftshift( kz.* ky  ./ (k2 + eps) );     

elseif mode == 1 % Discrete kernel

    
FOV = N.*spatial_res;
center = 1+N/2;
kx = 1:N(1);
ky = 1:N(2);
kz = 1:N(3);

kx = kx - center(1);
ky = ky - center(2);
kz = kz - center(3);

delta_kx = 1/FOV(1);
delta_ky = 1/FOV(2);
delta_kz = 1/FOV(3);


kx = kx * delta_kx;
ky = ky * delta_ky;
kz = kz * delta_kz;

kx = reshape(kx,[length(kx),1,1]);
ky = reshape(ky,[1,length(ky),1]);
kz = reshape(kz,[1,1,length(kz)]);

kx = repmat(kx,[1,N(2),N(3)]);
ky = repmat(ky,[N(1),1,N(3)]);
kz = repmat(kz,[N(1),N(2),1]);

k2 = -3+cos(2*pi*kx)+cos(2*pi*ky)+cos(2*pi*kz);
k2(k2==0) = eps;
kernel = 1/3 - (-1+cos(2*pi*kz)) ./ k2;


DC = (kx==0) & (ky==0) & (kz==0);
kernel(DC==1) = 0;
kernel = fftshift(kernel);
%     
% [ky,kx,kz] = meshgrid(-N(2)/2:N(2)/2-1, -N(1)/2:N(1)/2-1, -N(3)/2:N(3)/2-1);
% 
% kx = cos(2*pi*kx / spatial_res(1) / N(1));
% ky = cos(2*pi*ky / spatial_res(2) / N(2));
% kz = cos(2*pi*kz / spatial_res(3) / N(3));
% 
% k2 = kx + ky + kz;
% 
% R_tot = eye(3);
% 
% kernel = fftshift( 1/3 - (kx * R_tot(3,1) + ky * R_tot(3,2) + kz * R_tot(3,3)-1) ./ (k2 + eps-3) );    
    
elseif mode == 2 % Integrated Green function
    
FOV = N.*spatial_res;
center = 1+N/2;
    
kx = 1:N(1);
ky = 1:N(2);
kz = 1:N(3);

kx = kx - center(1);
ky = ky - center(2);
kz = kz - center(3);


% determine the step sizes delta_kx, delta_ky, delta_kz in k-space
delta_kx = N(1)/FOV(1);
delta_ky = N(2)/FOV(2);
delta_kz = N(3)/FOV(3);


kx = kx * delta_kx;
ky = ky * delta_ky;
kz = kz * delta_kz;

kx = reshape(kx,[length(kx),1,1]);
ky = reshape(ky,[1,length(ky),1]);
kz = reshape(kz,[1,1,length(kz)]);

kx = repmat(kx,[1,N(2),N(3)]);
ky = repmat(ky,[N(1),1,N(3)]);
kz = repmat(kz,[N(1),N(2),1]);

if nargin < 4
    theta = 0;
end

k2 = kx.^2+ky.^2+kz.^2;
k2(k2==0) = eps;

spatial_kernel = zeros(N);
for i=1:N(1)
    for j=1:N(2)
        for k=1:N(3)
            x(1) = kx(i,j,k)-0.5*delta_kx;
            x(2) = kx(i,j,k)+0.5*delta_kx;
            y(1) = ky(i,j,k)-0.5*delta_ky;
            y(2) = ky(i,j,k)+0.5*delta_ky;
            z(1) = kz(i,j,k)-0.5*delta_kz;
            z(2) = kz(i,j,k)+0.5*delta_kz;
            
            m = zeros([2 2 2]);
            for di = 1:2
                for dj = 1:2
                    for dk = 1:2
                        m(di,dj,dk) = (-1)^(di+dj+dk) * atan(x(di)*y(dj)/(z(dk)*sqrt(x(di)^2+y(dj)^2+z(dk)^2)));
                    end
                end
            end
            spatial_kernel(i,j,k) = -(1/(pi*4)) *sum(m(:));       
            
            
        end
    end
end

DC = (kx==0) & (ky==0) & (kz==0);
spatial_kernel(DC==1) = spatial_kernel(DC==1)+1/3;

kernel = fftn(fftshift(spatial_kernel));
    

end


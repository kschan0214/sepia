function [ background, kappa ] = fitmodels( img, mask, weight, phi_e )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
N = size(img);
se = strel('sphere',3);
mask3=imerode(mask,se);

    kappa = 0;%181; %9.35*phs_scale;
%     phi0 = 0;
% 
%     a1 = sum( weight(:).*mask3(:).*phi_e(:) );
%     a2 = sum( weight(:).*mask3(:).*phi_e(:).*phi_e(:) );
%     a3 = sum( weight(:).*mask3(:) );
%     f1 = sum( weight(:).*mask3(:).*phi_e(:).*z2(:) );
%     f2 = sum( weight(:).*mask3(:).*z2(:) );
%     det = a1*a1-a2*a3;
%     
%     phi0 = (a1*f1 - a2*f2)/det;
%     kappa = (-a3*f1 + a1*f2)/det;
    
    [ky,kx,kz] = meshgrid(-N(2)/2:N(2)/2-1, -N(1)/2:N(1)/2-1, -N(3)/2:N(3)/2-1);
    kx = (kx / max(abs(kx(:))));
    ky = (ky / max(abs(ky(:))));
    kz = (kz / max(abs(kz(:))));
    
A(1,1) = sum( weight(:).*mask3(:) );
A(2,2) = sum( weight(:).*mask3(:).*kx(:).*kx(:) );
A(3,3) = sum( weight(:).*mask3(:).*ky(:).*ky(:) );
A(4,4) = sum( weight(:).*mask3(:).*kz(:).*kz(:) );
A(5,5) = sum( weight(:).*mask3(:).*phi_e(:).*phi_e(:) );

A(1,2) = sum( weight(:).*mask3(:).*kx(:) );
A(2,1) = A(1,2);
A(1,3) = sum( weight(:).*mask3(:).*ky(:) );
A(3,1) = A(1,3);
A(1,4) = sum( weight(:).*mask3(:).*kz(:) );
A(4,1) = A(1,4);
A(1,5) = sum( weight(:).*mask3(:).*phi_e(:) );
A(5,1) = A(1,5);

A(2,3) = sum( weight(:).*mask3(:).*ky(:).*kx(:) );
A(3,2) = A(2,3);
A(2,4) = sum( weight(:).*mask3(:).*kz(:).*kx(:) );
A(4,2) = A(2,4);
A(2,5) = sum( weight(:).*mask3(:).*phi_e(:).*kx(:) );
A(5,2) = A(2,5);

A(3,4) = sum( weight(:).*mask3(:).*kz(:).*ky(:) );
A(4,3) = A(3,4);
A(3,5) = sum( weight(:).*mask3(:).*phi_e(:).*ky(:) );
A(5,3) = A(3,5);

A(4,5) = sum( weight(:).*mask3(:).*phi_e(:).*kz(:) );
A(5,4) = A(4,5);

F(1) = sum( weight(:).*mask3(:).*img(:) );
F(2) = sum( weight(:).*mask3(:).*img(:).*kx(:) );
F(3) = sum( weight(:).*mask3(:).*img(:).*ky(:) );
F(4) = sum( weight(:).*mask3(:).*img(:).*kz(:) );
F(5) = sum( weight(:).*mask3(:).*img(:).*phi_e(:) );

kappa = A\F';

    %phi0 = (f2-sum(weight(:).*mask(:).*phi_e(:)*kappa))/a3;
    
    background = kappa(5)*phi_e+kappa(4)*kz+kappa(3)*ky+kappa(2)*ky+kappa(1);

end


function [ kuw, luw, delta ] = iterlap_unwrap2( img, mask, weight, iters )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
N = size(img);

if (nargin<2)
    mask = ones(N);
    weight = mask;
    iters = 250;
end
if (nargin<3)
    weight = mask;
    iters = 250;
end
if (nargin<4)
    iters = 250;
end


delta = img;
luw = zeros(N);
residual = true;

% Iterative Laplacian unwrapping
for it = 1:iters
    %puL = unwrapLaplacian(mask.*delta,N,[1 1 1]);
   [kpuL, puL] = gduw( delta, mask, [1 1 1] );
    luw = luw+puL;
    delta = angle(exp(1i*(img-luw)));
    if max(mask(:).*delta(:))<pi/2 && min(mask(:).*delta(:)) > -pi/2
        %display(it);
        residual = false;
        break;
    end
end


linearBackground = zeros(N);

if residual == true
    % Remove linear gradients from residual
    [ky,kx,kz] = meshgrid(-N(2)/2:N(2)/2-1, -N(1)/2:N(1)/2-1, -N(3)/2:N(3)/2-1);
    kx = (kx / max(abs(kx(:))));
    ky = (ky / max(abs(ky(:))));
    kz = (kz / max(abs(kz(:))));

    se = strel('sphere',3);
    mask3=imerode(mask,se);
    
    kappaV = -20:0.2:20;
    mp = delta( N(1)/2,N(2)/2,N(3)/2);
    delta = delta-mp;
    for k = 1:201
        cost(k) = sum( mask3(:).*weight(:).*(1-cos( delta(:) - kappaV(k)*kz(:) )) );
    end
    [min_val,idx]=min(cost(:));
    kappaL(3) = kappaV(idx);
    delta = delta - kappaL(3)*kz;
    for k = 1:201
        cost(k) = sum( mask3(:).*weight(:).*(1-cos( delta(:) - kappaV(k)*ky(:) )) );
    end
    [min_val,idx]=min(cost(:));
    kappaL(2) = kappaV(idx);
    delta = delta - kappaL(2)*ky;
    for k = 1:201
        cost(k) = sum( mask3(:).*weight(:).*(1-cos( delta(:) - kappaV(k)*kx(:) )) );
    end
    [min_val,idx]=min(cost(:));
    kappaL(1) = kappaV(idx);
    delta = delta - kappaL(1)*kx;
    
    
    
    linearBackground = kappaL(3)*kz+ kappaL(2)*ky+ kappaL(1)*kx +mp;
    luw = luw + linearBackground;
    
    delta = angle(exp(1i*(img-luw)));
    
    kappaV = -2:0.02:2;
    mp = delta( N(1)/2,N(2)/2,N(3)/2);
    delta = delta-mp;
    for k = 1:201
        cost(k) = sum( mask3(:).*weight(:).*(1-cos( delta(:) - kappaV(k)*kz(:) )) );
    end
    [min_val,idx]=min(cost(:));
    kappaL(3) = kappaV(idx);
    delta = delta - kappaL(3)*kz;
    for k = 1:201
        cost(k) = sum( mask3(:).*weight(:).*(1-cos( delta(:) - kappaV(k)*ky(:) )) );
    end
    [min_val,idx]=min(cost(:));
    kappaL(2) = kappaV(idx);
    delta = delta - kappaL(2)*ky;
    for k = 1:201
        cost(k) = sum( mask3(:).*weight(:).*(1-cos( delta(:) - kappaV(k)*kx(:) )) );
    end
    [min_val,idx]=min(cost(:));
    kappaL(1) = kappaV(idx);
    delta = delta - kappaL(1)*kx;
    
    
    linearBackground = linearBackground+kappaL(3)*kz+ kappaL(2)*ky+ kappaL(1)*kx +mp;
    luw = luw +kappaL(3)*kz+ kappaL(2)*ky+ kappaL(1)*kx +mp;
    
    delta = angle(exp(1i*(img-luw)));
    
    kappaV = -1:0.01:1;
    mp = delta( N(1)/2,N(2)/2,N(3)/2);
    delta = delta-mp;
    for k = 1:201
        cost(k) = sum( mask3(:).*weight(:).*(1-cos( delta(:) - kappaV(k)*kz(:) )) );
    end
    [min_val,idx]=min(cost(:));
    kappaL(3) = kappaV(idx);
    delta = delta - kappaL(3)*kz;
    for k = 1:201
        cost(k) = sum( mask3(:).*weight(:).*(1-cos( delta(:) - kappaV(k)*ky(:) )) );
    end
    [min_val,idx]=min(cost(:));
    kappaL(2) = kappaV(idx);
    delta = delta - kappaL(2)*ky;
    for k = 1:201
        cost(k) = sum( mask3(:).*weight(:).*(1-cos( delta(:) - kappaV(k)*kx(:) )) );
    end
    [min_val,idx]=min(cost(:));
    kappaL(1) = kappaV(idx);
    delta = delta - kappaL(1)*kx;
    
    
    
    %display(kappaL);
    
    linearBackground = linearBackground+kappaL(3)*kz+ kappaL(2)*ky+ kappaL(1)*kx +mp;
    luw = luw +kappaL(3)*kz+ kappaL(2)*ky+ kappaL(1)*kx +mp;
    
    delta = angle(exp(1i*(img-luw)));
end

% Iterative unwrapping by 2pi steps, using the iterative laplacian as
% approximation
kuw =img;
for i = 1:150
    out_old = kuw;
    kuw = kuw + 2*pi*round( (luw - kuw)/(2*pi) );
    
    if sum(abs(out_old(:)-kuw(:))) < 1
        break;
    end
    
end



end


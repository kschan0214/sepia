function [SH3D,Residuals,b] = spherical_harmonic_shimming(img,mask,order)

dim = size(img);

%% for volume fitting
% find non-zero elements
Indices = find(mask(:,:,:));
% get Cartesian coordinate of the non-zero elements
[x1,y1,z1] = ind2sub(size(img(:,:,:)),Indices);
% solutions of the linear operation
R = img(Indices);
% Make sure have more solutions than unknowns
if length(Indices)>(3*order)^2
    
    % Spherical harmonic terms
    model = Create_model(x1,y1,z1,dim,order);
    % get coefficients of the SH terms
    b = pinv(model)*R;

    clear model
    clear R
    
    % compute the SH field with the full data
    Indices     = find(ones(dim));
    [x1,y1,z1]  = ind2sub(size(img(:,:,:)),Indices);
    model       = Create_model(x1,y1,z1,dim,order);
    Fit         = model*b;
    
    clear model
    
    % SH field in 3D
    SH3D = reshape(Fit,dim);
    
else

    SH3D = zeros(size(img(:,:,:)));
    
end

% non SH field
Residuals = img - SH3D;

end

%% SH model
function model = Create_model(x1,y1,z1,dim,order)
if order > 4
    order = 4;
end

N_terms = (order+1)^2;
model = double(zeros(length(x1),N_terms));

% zeroth-order SH component
if order >= 0
    model(:,1) = 1; % constant term
end
% first-order SH components
if order >= 1
    x = reshape(x1-dim(1)/2,length(x1),1);
    y = reshape(y1-dim(2)/2,length(x1),1);
    z = reshape(z1-dim(3)/2,length(x1),1);
    model(:,2) = x;     % x 
    model(:,3) = y;     % y 
    model(:,4) = z;     % z 
end
% second-order SH components
if order >= 2
    x2 = x.^2;
    y2 = y.^2;
    z2 = z.^2;
    model(:,5) = x2 - y2;               % x^2-y^2
    model(:,6) = x.*y;                	% xy 
    model(:,7) = 2*z2 - x2 - y2;        % 2z^2 - x^2 - y^2
    model(:,8) = y.*z;                 	% yz
    model(:,9) = x.*z;              	% xz
end
% third-order SH components
if order >= 3
    model(:,10) = y.*(3*x2 - y2);         	% y*(3x^2-y^2)
    model(:,11) = x.*y.*z;                	% xyz
    model(:,12) = y.*(4*z2 - x2 - y2);     	% y*(4z^2-x^2-y^2)
    model(:,13) = z.*(2*z2 - 3*x2 - 3*y2);  % z*(2z^2-3x^2-3y^2)
    model(:,14) = x.*(4*z2 - x2 - y2);   	% x(4z^2-x^2-y^2)
    model(:,15) = (x2 - y2).*z;          	% (x^2-y^2)z
    model(:,16) = (x2 - 3*y2).*x;          	% (x^2-3y^2)x
end
% fourth-order SH components
if order >= 4
    r   = sqrt(x2+y2+z2);
    r2  = r.^2;
    model(:,17) = x.*y.*(x2 - y2);                  % xy(x^2-y^2)
    model(:,18) = (3*x2 - y2).*y.*z;                % (3x^2-y^2)yz
    model(:,19) = x.*y.*(7*z2 - r2);                % xy(7z^2-3r2)
    model(:,20) = y.*z.*(7*z2 - 3*r2);              % yz(7z^2-3r^2)
    model(:,21) = 35*(z.^4) - 30*z2.*r2 + 3*(r.*4);	% 35z^4-30z^2r^2+3r^4
    model(:,22) = x.*z.*(7*z2 - 3*r2);              % xz(7z^2-3r^2)
    model(:,23) = (x2 - y2).*(7*z2 - r2);           % (x^2-y^2)(7z^2-r^2)
    model(:,24) = (x2 - 3*y2).*x.*z;                % (x^2-3y^2)xz
    model(:,25) = x2.*(x2 - 3*y2) - y2.*(3*x2 - y2); % x^2(x^2-3y^2) - y^2(3x^2-y^2)
end

end

% Generation of the Dipole Kernel
%   D = dipole_kernel(matrix_size, voxel_size, B0_dir, ...)
%
%   output
%   D - dipole kernel saved in Fourier space
% 
%   input
%   matrix_size - the size of the matrix
%   voxel_size - the size of a voxel
%   B0_dir - the direction of the B0 field
%   domain - 'imagespace' or 'kspace'
%       Fourier domain expression:
%       Salomir et al. Concepts in Magn Reson Part B 2003, 19(B):26-34
%       Image domain expression: 
%       Li et al. Magn Reson Med 2004, 51(5):1077-82
%
%
%   Created by Tian Liu in 2008
%   Modified by Tian Liu on 2011.02.01
%   Last modified by Tian Liu on 2013.07.22





function D=dipole_kernel(varargin)

[matrix_size voxel_size B0_dir domain] = parse_inputs(varargin{:});

if (B0_dir == 1)
    B0_dir = [1 0 0 ]';
elseif (B0_dir == 2)
    B0_dir = [0 1 0 ]';
elseif (B0_dir==3)
    B0_dir = [0 0 1]';
end

if strcmp(domain,'kspace')
    [Y,X,Z]=meshgrid(-matrix_size(2)/2:(matrix_size(2)/2-1),...
        -matrix_size(1)/2:(matrix_size(1)/2-1),...
        -matrix_size(3)/2:(matrix_size(3)/2-1));
    
    X = X/(matrix_size(1)*voxel_size(1));
    Y = Y/(matrix_size(2)*voxel_size(2));
    Z = Z/(matrix_size(3)*voxel_size(3));
    
    D = 1/3-  ( X*B0_dir(1) + Y*B0_dir(2) + Z*B0_dir(3) ).^2./(X.^2+Y.^2+Z.^2);
    D(isnan(D)) = 0;
    D = fftshift(D);
    
elseif strcmp(domain,'imagespace')
    [Y,X,Z]=meshgrid(-matrix_size(2)/2:(matrix_size(2)/2-1),...
        -matrix_size(1)/2:(matrix_size(1)/2-1),...
        -matrix_size(3)/2:(matrix_size(3)/2-1));
    
    X = X*voxel_size(1);
    Y = Y*voxel_size(2);
    Z = Z*voxel_size(3);
    
    d = (3*( X*B0_dir(1) + Y*B0_dir(2) + Z*B0_dir(3)).^2 - X.^2-Y.^2-Z.^2)./(4*pi*(X.^2+Y.^2+Z.^2).^2.5);

    d(isnan(d)) = 0;
    D = fftn(fftshift(d));
end


end


function [matrix_size voxel_size B0_dir domain] = parse_inputs(varargin)
    
    if size(varargin,2)<3
        error('At least matrix_size, voxel_size and B0_dir are required');
    end
    
    matrix_size = varargin{1};
    voxel_size = varargin{2};
    B0_dir = varargin{3};
    domain = 'kspace';
    
    if size(varargin,2)>3
        for k=4:size(varargin,2)
            if strcmpi(varargin{k},'imagespace')
                domain = 'imagespace';
            end
        end
    end
    
end

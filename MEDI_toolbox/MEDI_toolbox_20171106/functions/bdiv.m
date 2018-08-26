% Discrete Divergence Using Backward Difference
% with the Dirichlet Boundary Condition

% Created by Youngwook Kee (Oct 21 2015)
% Last modified date: Oct 24 2015

% References:
% [1] Chambolle, An Algorithm for Total Variation Minimization and
% Applications, JMIV 2004

function div = bdiv(Gx, voxel_size)

if (nargin < 2)
    voxel_size = [1 1 1];
end

% Gx = double(Gx);

Gx_x = Gx(:,:,:,1);
Gx_y = Gx(:,:,:,2);
Gx_z = Gx(:,:,:,3);

[Mx, My, Mz] = size(Gx_x);

Dx = [Gx_x(1:end-1,:,:); zeros(1,My,Mz)]...
    - [zeros(1,My,Mz); Gx_x(1:end-1,:,:)];

Dy = [Gx_y(:,1:end-1,:), zeros(Mx,1,Mz)]...
    - [zeros(Mx,1,Mz), Gx_y(:,1:end-1,:)];

Dz = cat(3, Gx_z(:,:,1:end-1), zeros(Mx,My,1))...
    - cat(3, zeros(Mx,My,1), Gx_z(:,:,1:end-1));

Dx = Dx/voxel_size(1);
Dy = Dy/voxel_size(2);
Dz = Dz/voxel_size(3);

div = -( Dx + Dy + Dz );

end

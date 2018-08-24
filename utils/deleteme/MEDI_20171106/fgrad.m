% Discrete Gradient Using Forward Differences
% with the Neuman Boundary Condition

% Created by Youngwook Kee (Oct 21 2015)
% Last modified date: Oct 24 2015

% References:
% [1] Chambolle, An Algorithm for Total Variation Minimization and
% Applications, JMIV 2004
% [2] Pock et al., Global Solutions of Variational Models with Convex
% Regularization, SIIMS 2010

function Gx = fgrad(chi, voxel_size)

if (nargin < 2)
    voxel_size = [1 1 1];
end

% chi = double(chi);

Dx = [chi(2:end,:,:); chi(end,:,:)] - chi;
Dy = [chi(:,2:end,:), chi(:,end,:)] - chi;
Dz = cat(3, chi(:,:,2:end), chi(:,:,end)) - chi;

Dx = Dx/voxel_size(1);
Dy = Dy/voxel_size(2);
Dz = Dz/voxel_size(3);

Gx = cat(4, Dx, Dy, Dz);

end

function [ phi_e ] = backmodel( airmask, body_rel_susc_brain, mode, padding, spatial_res )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
x_air_water = 9.395;
N = size(airmask);
if (nargin<2)
    body_rel_susc_brain = 0.0;
    mode = 2;
    padding = N;
    spatial_res = [1 1 1];
end
if (nargin<3)
    mode = 2;
    padding = N;
    spatial_res = [1 1 1];
end
if (nargin<4)
    padding = N;
    spatial_res = [1 1 1];
end
if (nargin<5)
    spatial_res = [1 1 1];
end
if length(padding) < 3
    padding(2) =  padding(1);
    padding(3) =  padding(1);
end

for z=1:N(3)
    slice = airmask(:,:,z);
    mval = max(slice(:));
    if mval == 1
        zadd = z;
        break
    end
end

for y=1:N(2)
    slice = squeeze(airmask(:,y,:));
    mval = max(slice(:));
    if mval == 1
        yadd = y;
        break
    end
end
for x=1:N(1)
    slice = squeeze(airmask(x,:,:));
    mval = max(slice(:));
    if mval == 1
        xadd = x;
        break
    end
end


pp = x_air_water*ones( N+2*padding );
pp((padding(1)+1):(padding(1)+N(1)),(padding(2)+1):(padding(2)+N(2)),(padding(3)+1):(padding(3)+N(3))) = x_air_water*(1-airmask);
pp((padding(1)+1+xadd):(padding(1)+N(1)-xadd),(padding(2)+1+yadd):(padding(2)+N(2)-yadd),1:(padding(3)+zadd)) = body_rel_susc_brain;
%imagesc3d2(pp, padding+N/2, 115, [90,90,-90], [-40,40], 0, 'air');

if mode == 2
    kernelg = dipole_kernel( N+2*padding,[1 1 1], 2 );
    pp = real(ifftn( kernelg.*fftn(pp) ));
else
    kernel = dipole_kernel( N+2*padding,[1 1 1], 0 );
    pp = real(ifftn( kernel.*fftn(pp) ));
end

phi_e = pp((padding(1)+1):(padding(1)+N(1)),(padding(2)+1):(padding(2)+N(2)),(padding(3)+1):(padding(3)+N(3)));

end


function erg = energy_ho_3dW2(voxel_size,iMag,kappa,psi,base_xy,base_z,p,cliques,disc_bar,disc_bar_z,th,quant)

% compute energy considering the magnitude weighting 

%energy_ho   Energy from kappa labeling and psi phase measurements.
%   erg = energy_ho(kappa,psi,base,p,cliques,disc_bar,p,th,quant) returns the energy of kappa labeling given the 
%   psi measurements image, the base ROI image (having ones in the region of interest (psi) and a passe-partout
%   made of zeros), the exponent p, the cliques matrix (each row indicating a displacement vector corresponding
%   to each clique), the disc_bar (complement to one of the quality maps), a threshold th defining a region for
%   which the potential (before a possible quantization) is quadratic, and quant which is a flag defining whether
%   the potential is or is not quantized.
%   (see J. Bioucas-Dias and G. Valadão, "Phase Unwrapping via Graph Cuts"
%   submitted to IEEE Transactions Image Processing, October, 2005).
%   SITE: www.lx.it.pt/~bioucas/ 

[m,n,h] = size(psi);
[cliquesm,cliquesn] = size(cliques); % Size of input cliques
maxdesl = max(max(abs(cliques))); % This is the maximum clique length used
% Here we put a passe-partout (constant length = maxdesl+1) in the images kappa and psi
base_kappa_xy = zeros(2*maxdesl+2+m,2*maxdesl+2+n,h);
for hh = 1:h
    base_kappa_xy(maxdesl+2:maxdesl+2+m-1,maxdesl+2:maxdesl+2+n-1,hh) = kappa(:,:,hh);
end
% To_do: the compute the base_kappa_z
base_kappa_z = zeros(m,n,2*maxdesl+2+h);

psi_base_xy      = zeros(2*maxdesl+2+m,2*maxdesl+2+n,h); 
iMag_base_xy = zeros(2*maxdesl+2+m,2*maxdesl+2+n,h);

for hh = 1:h
    psi_base_xy(maxdesl+2:maxdesl+2+m-1,maxdesl+2:maxdesl+2+n-1,hh) = psi(:,:,hh);
    iMag_base_xy(maxdesl+2:maxdesl+2+m-1,maxdesl+2:maxdesl+2+n-1,hh) = iMag(:,:,hh);        
end 
% To_do: the formulaiton of psi_base_z 

z = size(disc_bar,3);
base_disc_bar  = repmat(zeros(2*maxdesl+2+m,2*maxdesl+2+n),[1 1 z size(psi,3)]); 
base_disc_bar(maxdesl+2:maxdesl+2+m-1,maxdesl+2:maxdesl+2+n-1,:,:) = disc_bar;
base_disc_bar_z = disc_bar_z;

% compute the h and v energies; for zz = 1:height 
for zz = 1:size(psi,3)
    for t = 1:cliquesm
        % The allowed start and end pixels of the "interpixel" directed edge
        base_start(:,:,t,zz) = circshift(base_xy(:,:,t),[-cliques(t,1),-cliques(t,2)]).*base_xy(:,:,t);
        base_end(:,:,t,zz) = circshift(base_xy(:,:,t),[cliques(t,1),cliques(t,2)]).*base_xy(:,:,t);

        % By convention the difference images have the same size as the
        % original ones; the difference information is retrieved in the
        % pixel of the image that is subtracted (end of the diff vector)
        auxili(:,:,zz) = circshift(base_kappa_xy(:,:,zz),[cliques(t,1),cliques(t,2)]);
        t_dkappa(:,:,t,zz) = (base_kappa_xy(:,:,zz)-auxili(:,:,zz));   % the labels in big row minis one smaller row
        auxili2(:,:,zz) = circshift(psi_base_xy(:,:,zz),[cliques(t,1),cliques(t,2)]);
        
        dpsi(:,:,zz) = (auxili2(:,:,zz) - psi_base_xy(:,:,zz));
        
        auxili_iMag(:,:,zz) = circshift(iMag_base_xy(:,:,zz),[cliques(t,1),cliques(t,2)]);
        
        % note: we need to make the 257:260 to zero!
        d_iMag(:,:,t,zz)  = ((auxili_iMag(:,:,zz) + iMag_base_xy(:,:,zz))/2).*base_xy(:,:,zz).*circshift(base_xy(:,:,zz),[cliques(t,1),cliques(t,2)])...
                   .*base_disc_bar(:,:,t,zz);
                
        
        % Beyond base, we must multiply by
        % circshift(base,[cliques(t,1),cliques(t,2)]) in order to
        % account for frontier pixels that can't have links outside ROI
        a(:,:,t,zz) = (2*pi*t_dkappa(:,:,t,zz)-dpsi(:,:,zz)).*base_xy(:,:,zz).*circshift(base_xy(:,:,zz),[cliques(t,1),cliques(t,2)])...
                   .*base_disc_bar(:,:,t,zz);
    end
end

% To_do: to compute the energy of z directions 
base_kappa_z = kappa;
dkappa(:,:,1) = zeros(m,n);
dpsi_z(:,:,1) = zeros(m,n);
diMag_z(:,:,1) = zeros(m,n);
for zz = 2:size(psi,3)
    % the plane with higher z minus the plane with one smaller plane
    dkappa(:,:,zz) = base_kappa_z(:,:,zz) - base_kappa_z(:,:,zz-1);
    % the voxel potential in smaller z minus the ones with one bigger ones
    dpsi_z(:,:,zz) = psi(:,:,zz-1) - psi(:,:,zz);    
%     diMag_z(:,:,zz) = iMag(:,:,zz-1) - iMag(:,:,zz);  % wong in writte 
    diMag_z(:,:,zz) = (iMag(:,:,zz-1) + iMag(:,:,zz))/2;  
    b(:,:,zz) = (2*pi*dkappa(:,:,zz) - dpsi_z(:,:,zz)).*base_disc_bar_z(:,:,zz);
end

for zz = 1:size(psi,3)
    % when computing the energy, add voxel_size
    % assume that the voxel_size of x and y usually same. 
    erg0(zz) = sum(sum(sum(d_iMag(:,:,:,zz).*(clique_energy_ho(a(:,:,:,zz)/voxel_size(1),p,th,quant)))));
end
erg = sum(erg0);
% the energy of z direction
erg_b = sum(sum(sum(diMag_z.*(clique_energy_ho(b/voxel_size(3),p,th,quant)))));
erg = erg + erg_b;





% Projection onto Dipole Fields (PDF)
%   [RDF shim] = PDF(iFreq, N_std, Mask,matrix_size,voxel_size, B0_dir,tol)
% 
%   output
%   RDF - the relative difference field, or local field
%   shim (optional) - the cropped background dipole distribution
% 
%   input
%   iFreq - the unwrapped field map
%   N_std - the noise standard deviation on the field map. (1 over SNR for single echo)
%   Mask - a binary 3D matrix denoting the Region Of Interest
%   matrix_size - the size of the 3D matrix
%   voxel_size - the size of the voxel in mm
%   B0_dir - the direction of the B0 field
%   tol(optional) - tolerance level
%
%   When using the code, please cite 
%   T. Liu et al. NMR Biomed 2011;24(9):1129-36
%   de Rochefort et al. MRM 2010;63(1):194-206
%
%   Created by Tian Liu in 2009
%   Modified by Tian Liu on 2011.02.01
%   Last modified by Tian Liu on 2013.07.24


function [RDF shim] = PDF(iFreq, N_std, Mask, matrix_size, voxel_size, B0_dir, tol, n_CG, space, n_pad)

if (nargin<10)
    n_pad = 40;
end

if (nargin<9)
    space = 'imagespace';
end

if (nargin<8)
    n_CG = 30;
end

if (nargin<7)
    tol = 0.1;
end


% zero pad
matrix_size0 = matrix_size;
d1 = max(max(Mask,[],2),[],3);
d1first = find(d1,1,'first');
d1last = find(d1,1,'last');

d2 = max(max(Mask,[],1),[],3);
d2first = find(d2,1,'first');
d2last = find(d2,1,'last');
        
d3 = max(max(Mask,[],1),[],2);
d3first = find(d3,1,'first');
d3last = find(d3,1,'last');

if n_pad > 0
    matrix_size = [ floor((d1last - d1first+n_pad)/2)*2,...
                    floor((d2last - d2first+n_pad)/2)*2,...
                    floor((d3last - d3first+n_pad)/2)*2];
    iFreq=iFreq(d1first:d1last,d2first:d2last, d3first:d3last);
    N_std=N_std(d1first:d1last,d2first:d2last, d3first:d3last);
    Mask=Mask(d1first:d1last,d2first:d2last, d3first:d3last);
    padsize = [matrix_size(1)-size(iFreq,1) matrix_size(2)-size(iFreq,2) matrix_size(3)-size(iFreq,3)];
    iFreq = padarray(iFreq, padsize, 0,'post');
    N_std = padarray(N_std, padsize, 0,'post');
    Mask = padarray(Mask, padsize, 0,'post');
end


% generate the weighting
W = 1./N_std;
W(isinf(W)) =0;
W = W.*(Mask>0);
W_std = W;
W_var = W.^2;

%%%%%% start the PDF method %%%%%
if norm(B0_dir(:), 1) < 1.01
    D = dipole_kernel(matrix_size,voxel_size,B0_dir);
else % Oblique B0
    D = dipole_kernel(matrix_size,voxel_size,B0_dir,space);
end

% generating the RHS vector in Eq. 6 in the PDF paper
p_temp = real(ifftn(D.*fftn(W_var.*(iFreq) )));
b = p_temp( Mask(:) == 0);

% set erance level and maximum iteration allowed
E_noise_level = real(ifftn(D.*fftn(W_std.*ones(size(N_std)))));
itermax= n_CG;
A=@(xx)(dipole_term(W_var,D,Mask,xx) );
cg_tol = tol*norm(E_noise_level( Mask(:) == 0))/norm(b(:));

[x res num_iter] = cgsolve(A, b, cg_tol, itermax, 0);
fprintf('CG stops at: res %f, iter %d\n', res, num_iter);
 
xx = zeros(size(D));
xx(Mask(:) == 0) = x(1:end);
xx(Mask(:) > 0) = 0;

% background dipole field
p_dipole = real(ifftn(D.*fftn(xx)));


p_final = (iFreq-p_dipole).*Mask;  

% remove zero pad
if n_pad > 0
    RDF = zeros(matrix_size0);
    RDF(d1first:d1last,d2first:d2last, d3first:d3last) = ...
        p_final(1:d1last-d1first+1, 1:d2last-d2first+1, 1:d3last-d3first+1);

    shim = zeros(matrix_size0);
    shim(d1first:d1last,d2first:d2last, d3first:d3last) = ...
        xx(1:d1last-d1first+1, 1:d2last-d2first+1, 1:d3last-d3first+1);
else
    RDF = p_final;
    shim = xx;
end

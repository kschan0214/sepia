%% function chi = qsmSingleStepVSHARP(totalField,mask,matrixSize,voxelSize,varargin)
%
% Usage:
%   chi = qsmFANSI(totalField,mask,matrixSize,voxelSize,...
%           'tol',1e-2,'lambda',2.9e-2,'iteration',50,'magnitude',magn);
%
% Input
% --------------
%   totalField      : unwrapped total field
%   mask            : user-defined ROI mask
%   matrixSize      : image matrix size
%   voxelSize       : spatial resolution of image 
%   varargin        : flags with
%       'lambda'        -   user defined regularisation parameter for gradient L1 penalty
%       'fieldStrength'	-	magntic field strength of the scanner
%       'tol'           - tolerance for iteration
%       'iteration'     - maximum number of iterations
%       'vkernel'       - VSHARP kernel
%       (disabled)'magnitude'  	- magnitude images 
%       (disabled)'te'            - echo spacing
%
% Output
% --------------
%   chi             : QSM
%
% Description: single step QSM using VSHARP
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 15 July 2017
% Date last modified:
%
function chi = qsmSingleStepVSHARP(totalField,mask,matrixSize,voxelSize,varargin)
% parse input arguments
[lambda,magn,tol,maxiter,Kernel_Sizes] = parse_vararginSSQSM(varargin);
% [B0,TE,lambda,magn,tol,maxiter,Kernel_Sizes] = parse_vararginSSQSM(varargin);
% gyro = 2*pi*42.58;

% total field in k-space
kTotalField = fftn(totalField);

% create dipole kernel
D = DipoleKernel(matrixSize,voxelSize);
% gradient masks from magnitude image using k-space gradients
[fdx,fdy,fdz,E2] = GradientOperatorKspace(matrixSize);
cfdx = conj(fdx);       cfdy = conj(fdy);       cfdz = conj(fdz);

%% V-Sharp filtering for background removal
DiffMask = zeros([matrixSize, length(Kernel_Sizes)]);
Mask_Sharp = zeros([matrixSize, length(Kernel_Sizes)]);
Del_Sharp = zeros([matrixSize, length(Kernel_Sizes)]);
phs = 0;
for k = 1:length(Kernel_Sizes)
    Kernel_Size = Kernel_Sizes(k);

    ksize = [Kernel_Size, Kernel_Size, Kernel_Size];                % Sharp kernel size

    khsize = (ksize-1)/2;
    [a,b,c] = meshgrid(-khsize(2):khsize(2), -khsize(1):khsize(1), -khsize(3):khsize(3));

    kernel = (a.^2 / khsize(1)^2 + b.^2 / khsize(2)^2 + c.^2 / khsize(3)^2 ) <= 1;
    kernel = -kernel / sum(kernel(:));
    kernel(khsize(1)+1,khsize(2)+1,khsize(3)+1) = 1 + kernel(khsize(1)+1,khsize(2)+1,khsize(3)+1);

    Kernel = zeros(matrixSize);
    Kernel( 1+matrixSize(1)/2 - khsize(1) : 1+matrixSize(1)/2 + khsize(1), ...
            1+matrixSize(2)/2 - khsize(2) : 1+matrixSize(2)/2 + khsize(2), ...
            1+matrixSize(3)/2 - khsize(3) : 1+matrixSize(3)/2 + khsize(3) ) = kernel;

    del_sharp = fftn(fftshift(Kernel));    

    % erode mask to remove convolution artifacts
    erode_size = ksize + 1;

    msk_sharp = imerode(mask, strel('line', erode_size(1), 0));
    msk_sharp = imerode(msk_sharp, strel('line', erode_size(2), 90));
    msk_sharp = permute(msk_sharp, [1,3,2]);
    msk_sharp = imerode(msk_sharp, strel('line', erode_size(3), 0));
    msk_sharp = permute(msk_sharp, [1,3,2]);

    Mask_Sharp(:,:,:,k) = msk_sharp; 
    Del_Sharp(:,:,:,k) = del_sharp; 
    
    if k == 1
        DiffMask(:,:,:,1) = Mask_Sharp(:,:,:,1);
    else
        DiffMask(:,:,:,k) = Mask_Sharp(:,:,:,k) - Mask_Sharp(:,:,:,k-1);
    end
    
%     phs = phs + DiffMask(:,:,:,k) .* ifftn(Del_Sharp(:,:,:,k) .* kTotalField) / (TE * B0 * gyro);
        phs = phs + DiffMask(:,:,:,k) .* ifftn(Del_Sharp(:,:,:,k) .* kTotalField);
end
mask_sharp = Mask_Sharp(:,:,:,end);     % largest mask

%% gradient masks from magnitude image using k-space gradients
magni = magn .* mask_sharp;

Magn = fftn(magni / max(magni(:)));
magn_grad = cat(4, ifftn(Magn.*fdx), ifftn(Magn.*fdy), ifftn(Magn.*fdz));

magn_weight = zeros(size(magn_grad));
for s = 1:size(magn_grad,4)
    magn_use = abs(magn_grad(:,:,:,s));
    
    magn_order = sort(magn_use(mask_sharp==1), 'descend');
    magn_threshold = magn_order( round(length(magn_order) * .15) );
    magn_weight(:,:,:,s) = magn_use <= magn_threshold;
end

%% Single-step QSM with V-Sharp
Rhs_pcg = 0;
for k = 1:size(Del_Sharp,4)
    Rhs_pcg = Rhs_pcg + conj(Del_Sharp(:,:,:,k)) .* fftn(DiffMask(:,:,:,k) .* phs);
end
Rhs_pcg = conj(D) .* Rhs_pcg;       % right hand side

B_inv = 1 ./ (eps + abs(Del_Sharp(:,:,:,1) .* D).^2 + lambda*E2);        % preconditioner
precond_inv = @(x, B_inv) B_inv(:).*x;

x0 = B_inv .* Rhs_pcg;              % initial guess

tic
    [F_chi, flag, pcg_res, pcg_iter] = pcg(@(x)SSQSM_vsharp(x, D, Del_Sharp, conj(Del_Sharp), DiffMask, lambda, fdx, fdy, fdz, cfdx, cfdy, cfdz, magn_weight), ...
        Rhs_pcg(:), tol, maxiter, @(x) precond_inv(x, B_inv), [], x0(:));
toc

disp(['PCG iter: ', num2str(pcg_iter), '   PCG residual: ', num2str(pcg_res)])

chi = real(ifftn(reshape(F_chi, matrixSize))) .* mask_sharp;
end

%% Parse input arguments
% function [lambda,magn,tol,maxiter,Kernel_Sizes]=parse_vararginSSQSM(arg)
% % function [B0,TE,lambda,magn,tol,maxiter,Kernel_Sizes]=parse_vararginSSQSM(arg)
% % B0 = 3;
% % TE = 1;             %second
% lambda = 2.9e-2;
% magn = [];
% maxiter = 30;
% tol = 1e-2;
% Kernel_Sizes = 11:-2:3;
% 
% if ~isempty(arg)
%     for kvar = 1:length(arg)
% %         if strcmpi(arg{kvar},'fieldStrength')
% %             B0 = arg{kvar+1};
% %         end
% %         if strcmpi(arg{kvar},'te')
% %             TE = arg{kvar+1};
% %         end
%         if strcmpi(arg{kvar},'tol')
%             tol = arg{kvar+1};
%         end
%         if strcmpi(arg{kvar},'iteration')
%             maxiter = arg{kvar+1};
%         end
%         if strcmpi(arg{kvar},'magnitude')
%             magn = arg{kvar+1};
%         end
%         if strcmpi(arg{kvar},'lambda')
%             lambda = arg{kvar+1};
%         end
%         if strcmpi(arg{kvar},'vkernel')
%             if ~isempty(arg{kvar+1})
%                 Kernel_Sizes = arg{kvar+1};
%             end
%         end
%     end
% end
% end
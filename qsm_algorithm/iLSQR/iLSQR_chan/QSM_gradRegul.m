function [Ax] = QSM_gradRegul(X_vector,params_in,tflag)
% X_vector = vectorized Matrix of x (susceptibility in the case where the direct )
%
% params_in has the following fields
% kernel - k-space representation of Dipole Kernel
% msk = mask represent the regions where
% dims = 3d-dimension of the field map
% Ex,Ey,Ez - 3D K-space partial derivative operators 
% lambda = regularisation parameter of weighting matrices
%
% transp_flag = parameter for lsqr solver
% ------------------------------------------------------------------------

Ex = params_in.Ex;
Ey = params_in.Ey;
Ez = params_in.Ez;

if strcmp(tflag,'transp')      % y = A'*x
    % KC: Input x is Mx1 from 'notransp' A*x, where A is MxN and x is Nx1. 
    %     A' has dimension NxM such that A'*x is Nx1; A' is conjugate 
    %     transpose of A. The A.' will change the directions of submatrices 
    %     from vertical to horizontal, i.e. A = [DK;Ex;Ey;Ez]; 
    %     A.' = [DK,Ex,Ey,Ez], since they are all diagonal matrices.
    %     A' can be achieved by conj([DK,Ex,Ey,Ez]) 
    %     The whole operation can be broken down into 4 submatrices
    %     operations:  
    
    %   reshape the X_vector
    % KC: X_vector consists results of 4 submatrices
    X_matrix=reshape(X_vector,[params_in.dims 4]);
    
    %   apply the weighting mask the four submatrices
    X = bsxfun(@times,params_in.msk, X_matrix);
    
%     if true
%         ATx = zeros([params_in.dims 4], 'gpuArray');
%     else
        ATx = zeros([params_in.dims 4]);
%     end
    % --- calculation of FT^{-1}( conj(kernel) *FT(X))
    % KC: work on QSM submatrix problem
    ATx(:,:,:,1) = ifftn( conj(params_in.kernel) .* fftn(X(:,:,:,1)) );
        
    % --- apply partial derivative along x
    % KC: in r-space
    ATx(:,:,:,2) = params_in.lambda*ifftn( conj(Ex).*fftn(X(:,:,:,2)) );
   
    % --- apply partial derivative along y
    % KC: in r-space
    ATx(:,:,:,3) = params_in.lambda*ifftn( conj(Ey).*fftn(X(:,:,:,3)) );
    
    % --- apply partial derivative along z
    % KC: in r-space
    ATx(:,:,:,4) = params_in.lambda*ifftn( conj(Ez).*fftn(X(:,:,:,4)) );
        
    % --- combine the various submatrices
    % KC: combination of the four sub-operations is done by element-wise
    %     addition such that the result is equivalent to A'*x
    Ax = sum(ATx,4);
        
    % --- our problem is bound to be real in real space
    Ax = real(Ax(:));
    
elseif strcmp(tflag,'notransp') % y = A*x
    
    % KC: susceptibility
    X_matrix_nt=reshape(X_vector,params_in.dims);

    %   calculation of Msk*FT^{-1}(C*FT(X))
    % KC: susceptibilty in k-space
    X_matrix_nt=fftn(X_matrix_nt);
    % KC: masked, estimated local field
    Ax_0 = params_in.msk .* ifftn( params_in.kernel .* X_matrix_nt );

    % --- application of the gradient operator and regularization parameter
    % KC: masked, gradient in r-space with regularisation
    Ax_1 = params_in.lambda*params_in.msk.*ifftn( Ex.*X_matrix_nt );
    Ax_2 = params_in.lambda*params_in.msk.*ifftn( Ey.*X_matrix_nt );
    Ax_3 = params_in.lambda*params_in.msk.*ifftn( Ez.*X_matrix_nt );
    
    % --- combine the various submatrices
    % KC: [Eq. 9]
    Ax = cat( 4 , Ax_0 , Ax_1 , Ax_2 , Ax_3 );
    % --- our problem is bound to be real in real space
    % KC: Ax has the same number of elements as b with order:
    %     local field>x-gradient>y-gradient>z-gradient
    Ax = real(Ax(:));
    
    fprintf('+');
end
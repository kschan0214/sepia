%% function [chi, lambdaOptimal] = qsmClosedFormL2(localField,mask,matrixSize,voxelSize,varargin) 
%
% Description: compute QSM based on closed-form solution
% Ref        : Bilgic et al. JMRI 40:181-191(2014)
%            : Bilgic et al. MRM 72:1444-1459(2014)
%
% Input
% -----
%   localField      : local field perturbatios
%   mask            : user-defined mask
%   matrixSize      : image matrix size
%   voxelSize       : spatial resolution of image 
%   varargin        : flags with
%       'lambda'    -   user define regularisation parameter
%       'optimise'  -	self-define regularisation based on curvature of 
%                       L-curve (ref.2)
% 
% Ouput
% -----
%   chi             : QSM
%   lamdaOptimal    : optimal regularisation value based on L-curve
% 
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 24 March 2017
% Date last modified: 6 September 2017
%
function [chi, lambdaOptimal] = qsmClosedFormL2(localField,mask,matrixSize,voxelSize,varargin)
DEBUG=false;
%% Parsing varargin
[lambda, optimise] = parse_varargin_CFL2norm(varargin);

% display message
if ~optimise
    disp('Self-defined regularisation parameter? No');
    fprintf('Regularisation parameter: %f \n',lambda);
else
     disp('Self-defined regularisation parameter? Yes');   
end

% dipole kernel
kernel = DipoleKernel(matrixSize,voxelSize);

%% core
% defining gradient operators in k-space
[~,~,~,EtE] = GradientOperatorKspace(matrixSize);

DtD = abs(kernel).^2;

% Avoid comput FFT multiple times for self-defined method
kLocalField = fftn(localField);

%% closed-form solution for QSM
switch optimise 
    case false
        % remove the conjugate operator has no effect in value since kernel is
        % real but this makes the code consistence with Eq.8 
        % chi = real( ifftn(conj(kernel) .* fftn(localField) ./ (DtD + lambda^2 * EtE))) .* mask;
        chi_cplx = ifftn(kernel .* kLocalField ./ (DtD + lambda^2 * EtE));
        lambdaOptimal = [];

    case true
        %% KC: self-defined lambda based on MRM 72:1444-1459 (2014)
        % multigrid approach, 3 iterations
        % start with coarse grid
        lambdaMin = 1e-3; lambdaMax = 1e-1;
        % KC: create linear-spaced lambda condidates
        lambdaCandidate = linspace(lambdaMin,lambdaMax,25);
        % KC: modify to unequal spacing
        lambdaCandidate = cumsum(lambdaCandidate);

        % KC: compute the residual norm for each lambda being used
        %     chi results won't be saved at this stage to avoid memory usage
        normDataFidelity = zeros(1,length(lambdaCandidate));
        normRegularisation = zeros(1,length(lambdaCandidate));
        for klambda = 1:length(lambdaCandidate)
            chi_temp = ifftn(kernel .* kLocalField ./ (DtD + lambdaCandidate(klambda)^2 * EtE));
            % KC: norm of Residual of data;
            [normDataFidelity(klambda),~] = ResidualGivenQSMLocalFieldDipolekernel(chi_temp, localField, kernel);
            % KC: norm of Gradients
            [normRegularisation(klambda), ~] = ResidualGradientGivenQSM(chi_temp,matrixSize);
        end

        % KC: finding optimal lambda based on curvature of L-curve  
        rho = log(normDataFidelity.^2);
        omega = log(normRegularisation.^2);
        % KC: First derivatives w.r.t. lambda
        drho = myDerivative(rho,lambdaCandidate);
        domega = myDerivative(omega,lambdaCandidate);
        % KC: Second derivatives w.r.t. lambda
        d2rho = myDerivative(drho,lambdaCandidate);
        d2omega = myDerivative(domega,lambdaCandidate);
        % KC: Eq. 19 of Bilgic et al. MRM 72:1444-1459 (2014)
        curvature = (d2rho.*domega - d2omega.*drho)./ (drho.^2 + domega.^2).^1.5;
        % use interpolation to have finer resolution of the curvature
        lambdaCandidateInterp = linspace(min(lambdaCandidate),max(lambdaCandidate),1000);
        curvatureInterp = interp1(lambdaCandidate,curvature,lambdaCandidateInterp,'spline');
        % KC: optimal when the curvature is maximum
        [~, I] = sort(curvatureInterp,'descend');
        lambdaOptimal = lambdaCandidateInterp(I(1));
        
        if DEBUG
            % KC: display L-curve and L-curve curvature and optimal lambda
            figure;plot(normDataFidelity,normRegularisation,'bx-');title('L-curve');
            xlabel('Data fidelity');ylabel('chi maps Gradient');
            figure;plot(lambdaCandidateInterp,curvatureInterp,'ro-');title('Curvature of L-curve');
            xlabel('Regularisation parameter');ylabel('Curvature');
            drawnow
        end
        display(['Curvature of L-curve is maximal when lambda = ' num2str(lambdaOptimal)]);
        
        % KC: apply closed-form solution to optimal lambda
        chi_cplx = ifftn(kernel .* kLocalField ./ (DtD + lambdaOptimal^2 * EtE));

end

% KC: return real part only
chi = real(chi_cplx.* mask);

end

%% KC: return derivative of x with respect to y using gradient function
% Use of bulit-in function diff will lose the no. of elements 
function res = myDerivative(x,y)
    res = gradient(x)./gradient(y);
end

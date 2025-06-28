%% [r2s,t2s,m0] = R2star_NLLS(img,te,mask,isParallel,NUM_MAGN)
%
% Input
% --------------
% img           : multiecho images, time in last dimension
% te            : echo times
% mask          : signal mask
% isParallel    : parallel computing (i.e. parfor)
% NUM_MAGN      : no. of phase-corrupted echoes
%
% Output
% --------------
% r2s           : R2* map
% t2s           : T2* map
% m0            : M0, T1-weighted image
%
% Description: R2* fitting by nonlinear least-squares(NLLS), parallel
% computing compatible
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: July 17, 2016
% Date last modified: 24 Sep 2017
%
function [r2s,t2s,s0] = R2star_NLLS(img,te,mask,isParallel,NUM_MAGN)

% set range of R2* and T2*
minT2s = min(te)/20;
maxT2s = max(te)*20;
ranget2s = [minT2s, maxT2s];
ranger2s = [1/maxT2s, 1/minT2s];

% if no mask input use default mask to speed up the fitting
if isempty(mask)
    mask = max(img,[],4)/(max(img(:))) >0.05;
else
    mask = mask > 0;
end

if nargin<4
    isParallel = false;
end

% check if the input is real then switch to magnitude fitting
if nargin<5 || isreal(img)
    NUM_MAGN = length(te);
end

% define number of phase corrupted echoes
if NUM_MAGN==length(te)
	disp('Magnitude fitting');
elseif NUM_MAGN==0
	disp('Complex fitting');
else
	disp('Mixed fitting');
	disp(['No. of phase corrupted echoes: ' num2str(NUM_MAGN)]);
end

% get dimension of data
matrixSize = size(img);

%% main
% estimate initial guesses with fast method
[r2s0,~,m00] = R2star_trapezoidal(img,te);

% convert initial guesses to 1D array
r2s0 = r2s0(:);
m00 = m00(:);
mask = mask(:);
img = reshape(img,[numel(img)/matrixSize(end), matrixSize(end)]);

% fitting algorithm parameters(boundaries)
lb = double([ranger2s(1), abs(min(img(:)))]);
ub = double([ranger2s(2), 2*max(img(:))*exp(ranger2s(2)*min(te))]);
options = optimoptions(@lsqnonlin,'Display','off','MaxIter',100);

r2s = zeros(size(r2s0));
s0 = zeros(size(m00));
if isParallel
    parfor k = 1:length(r2s0)
        if mask(k)
            % initial guesses
            c0 = double([r2s0(k), m00(k)]);

            % fitting core
            [res, ~] = lsqnonlin(@(x)fitError_R2s(x,squeeze(img(k,:)),te,NUM_MAGN), c0,lb,ub,options);

            % get results
            r2s(k) = res(1);
            s0(k) = res(2);
        else
            r2s(k) = 0;
            s0(k) = 0;
        end
    end
%     catch
else
    for k = 1:length(r2s0)
        if mask(k)
            % initial guesses
            c0 = double([r2s0(k), m00(k)]);

            % fitting core
            [res, ~] = lsqnonlin(@(x)fitError_R2s(x,squeeze(img(k,:)),te,NUM_MAGN), c0,lb,ub,options);

            % get results
            r2s(k) = res(1);
            s0(k) = res(2);
        else
            r2s(k) = 0;
            s0(k) = 0;
        end
    end
end

%% reshape 1D results to images
r2s = reshape(r2s,matrixSize(1:end-1));
s0 = reshape(s0,matrixSize(1:end-1));
t2s = 1./r2s;
t2s = SetImgRange(t2s,ranget2s);

end

%% Compute residual
function [fiter] = fitError_R2s(x,sMeas,te,NUM_MAGN)

r2s = x(1);
m0 = x(2);
sFit = Signal_mGRE(m0,r2s,te,'r2s');

% compute fitting residual
fiter = sepia_computeFiter(sFit(:),sMeas(:),NUM_MAGN);

end


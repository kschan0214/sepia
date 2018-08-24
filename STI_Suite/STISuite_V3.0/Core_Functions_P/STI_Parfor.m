function [SusceptibilityTensor]=STI_Parfor(Thetai,HVector,B0,TE_ms,parallel_Flag)
%% susceptibility tensor calculation
% Inputs: 
% Phase4D: 4D dataset the 4th dimention is orientation. e.g. 
%          size(Phase4D)=[128,128,128,n];                   % n directions
% H_Matrix: the H matrix
%          e.g. H= [  0.0072   -0.8902   -0.4556            % direction 1
%                     0.1121   -0.2320   -0.9662            % direction 2
%                     ........
%                    -0.1384    0.7923    0.5942];          % direction n
% B0: B0, e.g. B0=3;
% TE_ms: TE in the unit of ms, e.g TE=30.
% parallel_Flag: parallel computing flag. recommend 'off'.
% for example
%[SusceptibilityTensor]=STI_Parfor(Phase4D,H_Matrix,B0,TE_ms,parallel_Flag);
% Chunlei Liu, PHD
% UC Berkeley.

poolnum=1;
switch parallel_Flag
    case 'on'
        matlabpool local 3
        poolnum=3;
end

tic
xres=size(Thetai,1);
yres=size(Thetai,2);
zres=size(Thetai,3);
Thetai=padarray(Thetai,[0 xres-yres xres-zres 0]/2);

SZ1=xres^3;
N_direction=size(Thetai,4);
thetai0=zeros(SZ1,N_direction);
for i=1:N_direction
    Xk0=fftnc(Thetai(:,:,:,i));
    thetai0(:,i)=Xk0(:);
end
disp('data FFT_ed')
clear Xk0
H=HVector./repmat(sqrt(sum(HVector.^2,2)),[1,3]);
A1 =[H(:,1).^2  2*H(:,1).*H(:,2)  2*H(:,1).*H(:,3) H(:,2).^2 2*H(:,2).*H(:,3) H(:,3).^2]/3;
SS=[xres xres xres];
[ry,rx,rz] = meshgrid(-SS(2)/2:SS(2)/2-1,-SS(1)/2:SS(1)/2-1,-SS(3)/2:SS(3)/2-1);
ki=[rx(:) ry(:) rz(:)];
clear rx ry rz
disp('ki constructed')
%% Variable distributing
SamplingVector1=(round((SZ1/poolnum)*(1:poolnum)))';
SamplingVector0=[1; SamplingVector1(1:end-1)+1];

for i=1:poolnum
kipar{i}=ki(SamplingVector0(i):SamplingVector1(i),:);
thetai0par{i}=thetai0(SamplingVector0(i):SamplingVector1(i),:);
end

clear ki thetai0
disp('Parallel variables created')

%% Parfor loop
parfor loopvar = 1:poolnum
    warning off all
    xkipar{loopvar}=myparforfun(kipar{loopvar},A1,H,thetai0par{loopvar},loopvar);
end
clear kipar

xki=[];
for i=1:poolnum
    xki=[xki; xkipar{i}];
end
size(xki)
%% Organize
xki=reshape(xki,[xres,xres,xres,6]);

xki(isnan(xki))=0;
for i=1:6
    X_tensor(:,:,:,i)=ifftnc(xki(:,:,:,i));
end
disp('reconstruction done')

SF=ScalingFactor(B0,TE_ms);
SusceptibilityTensor=X_tensor*SF.X;

switch parallel_Flag
    case 'on'
     matlabpool close
end
function xki=myparforfun(ki,A1,H,thetai0,loopvar)
SZ1=size(ki,1);
tic
xki=zeros(SZ1,6);
for i=1:SZ1
    k=ki(i,:);
    A2 =[k(1)*H(:,1) k(1)*H(:,2)+k(2)*H(:,1) k(1)*H(:,3)+k(3)*H(:,1) k(2)*H(:,2) k(2)*H(:,3)+k(3)*H(:,2) k(3)*H(:,3)];
    A=A1-repmat(H*k',[1,6]).*A2/sum(k.^2);
    thetai= thetai0(i,:)';
    xki(i,:)=(A\thetai)';
    if ~mod(i,round(SZ1/20))
        disp([ num2str(loopvar) ': ' num2str(round(i/SZ1*100)) '% finished;   Remaining: ' num2str(round((1-i/SZ1)/(i/SZ1)*toc)) ' sec' ])
    end
end





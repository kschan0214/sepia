function [FIT3D,Residuals,b] = Shim2ndOrder(img,Mask,order)

% dim=Volumenii.hdr.dime.dim(2:4);
dim = size(img);

%% for volume fitting


Indices=find(Mask(:,:,:));
[x1,y1,z1]=ind2sub(size(img(:,:,:)),Indices);
R=img(Indices);
% keyboard
if length(Indices)>(3*order)^2
    
    model=Create_model(x1,y1,z1,dim,order);
%     keyboard
    %     b=model\R;
    b=pinv(model)*R;
%     temp=R-model*b;
    clear model
    clear R
%     [FIT3D,Residuals,b1,model1]=PolyFitNifti(Volumenii,Mask,order);
%     Residuals=Volumenii;
%     Residuals.img=Residuals.img*0;
%     Residuals.img(Indices)=temp;
%     subplot(121);Orthoview (Residuals.img)
%     subplot(122);Orthoview (Volumenii.img.*Mask.img)
    Indices = find(ones(dim));
    [x1,y1,z1]=ind2sub(size(img(:,:,:)),Indices);
    model=Create_model(x1,y1,z1,dim,order);
%     for k=1:10,subplot(5,2,k),plot(model(:,k)),end
    Fit=model*b;
    clear model
%     FIT3D=Volumenii;
    FIT3D=reshape(Fit,dim);
%     Residuals=FIT3D;
    Residuals=(img-FIT3D);
else
%     FIT3D=Volumenii;
    FIT3D=zeros(size(img(:,:,:)));
%     Residuals=FIT3D;
    Residuals=(img-FIT3D).*Mask;
end;

% %     figure
% %     keyboard
% %
% %     imagesc(FIT3D());colorbar
function model = Create_model(x1,y1,z1,dim,order);
Nsize=[1 4 10 20 35];
N=Nsize(order+1);
model=double(zeros(length(x1),N));
%zeroth order
if order>=0
    model(1:length(x1),1)=1;
end
%first order
if order>=1
    model(:,2)=reshape(x1-dim(1)/2,length(x1),1);%x -siemens
    model(:,3)=reshape(y1-dim(2)/2,length(x1),1);%y -siemens
    model(:,4)=reshape(z1-dim(3)/2,length(x1),1);%z -siemens
end
%second order
if order>=2

%    model(:,5) =model(:,2).*model(:,2);%x^2y^0z^0
%    model(:,6) =model(:,2).*model(:,3);%x^1y^1z^0 -siemens xy
%    model(:,7) =model(:,3).*model(:,3);%x^0y^2z^0
%    model(:,8) =model(:,3).*model(:,4);%x^0y^1z^1 -siemens yz
%    model(:,9) =model(:,4).*model(:,4);%x^0y^0z^2 -siemens 
%    model(:,10)=model(:,4).*model(:,2);%x^1y^0z^1 -siemens xz
    model(:,5) = model(:,2).*model(:,2) ...
				- model(:,3).*model(:,3);% siemens x^2-y^2
    model(:,6) = model(:,2).*model(:,3);%x^1y^1z^0 -siemens xy
    model(:,7) = 2 * model(:,4).*model(:,4) ...
				- model(:,3).*model(:,3)...
				- model(:,2).*model(:,2); % 2 z^2 - x^2 - y^2
    model(:,8) =model(:,3).*model(:,4);%x^0y^1z^1 -siemens yz
    model(:,9)=model(:,4).*model(:,2);%x^1y^0z^1 -siemens xz
end
%%third order
%if order>=3
%    model(:,11) =model(:,5).*model(:,2); %x^3y^0z^0
%    model(:,14) =model(:,7).*model(:,3); %x^0y^3z^0
%    model(:,18) =model(:,9).*model(:,4); %x^0y^0z^3
%    model(:,15) =model(:,8).*model(:,2); %x^1y^1z^1
%    model(:,12) =model(:,6).*model(:,2); %x^2y^1z^0
%    model(:,13) =model(:,7).*model(:,2); %x^1y^2z^0
%    model(:,16) =model(:,9).*model(:,2); %x^1y^0z^2
%    model(:,17) =model(:,10).*model(:,2);%x^2y^0z^1
%    model(:,19) =model(:,9).*model(:,3); %x^0y^1z^2
%    model(:,20) =model(:,8).*model(:,3); %x^0y^2z^1
%end
%if order>=4
%    model(:,21) =model(:,11).*model(:,2); %x^4y^0z^0
%    model(:,22) =model(:,14).*model(:,3); %x^0y^4z^0
%    model(:,23) =model(:,18).*model(:,4); %x^0y^0z^4
%    model(:,24) =model(:,11).*model(:,3); %x^3y^1z^0
%    model(:,25) =model(:,11).*model(:,4); %x^3y^0z^1
%    model(:,26) =model(:,5).*model(:,7); %x^2y^2z^0
%    model(:,27) =model(:,5).*model(:,9); %x^2y^0z^2
%    model(:,28) =model(:,17).*model(:,3);%x^2y^1z^1
%    model(:,29) =model(:,14).*model(:,4); %x^0y^3z^1
%    model(:,30) =model(:,14).*model(:,2); %x^1y^3z^0
%    model(:,31) =model(:,13).*model(:,4); %x^1y^2z^1
%    model(:,32) =model(:,20).*model(:,4); %x^0y^2z^2
%    model(:,33) =model(:,18).*model(:,3); %x^0y^1z^3
%    model(:,34) =model(:,18).*model(:,2); %x^1y^0z^3
%    model(:,35) =model(:,15).*model(:,4); %x^1y^1z^2
%end

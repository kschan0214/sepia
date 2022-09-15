function [FIT3D,Residuals,b]=PolyFit_weights(Data3D,Mask,order,weights);

dim=size(Data3D);


%% for volume fitting


Indices = find(Mask(:,:,:));
[x1,y1,z1]=ind2sub(size(Data3D(:,:,:)),Indices);
%     x1 = cast(x1,'like',Data3D);
%     y1 = cast(y1,'like',Data3D);
%     z1 = cast(z1,'like',Data3D);
if isempty(weights)
    w = ones([length(Indices),1]);
else
    w = weights(Indices);
end

R=Data3D(Indices);
model=Create_model(x1,y1,z1,dim,order);

if length(Indices)>size(model,2)
    %         b=model\R;
    %     b=pinv(model)*R;
    % the fit is weighted by the weigthing matrix
    b=pinv(bsxfun(@times,model,w))*(R.*w);
    Fit=model*b;
    clear model
    Indices = find(ones(dim));
    [x1,y1,z1]=ind2sub(dim,Indices);
    model=Create_model(x1,y1,z1,dim,order);
    Fit=model*b;
    FIT3D=reshape(Fit,dim);
    clear Fit
    Residuals=FIT3D;
    Residuals=(Data3D-FIT3D).*Mask;
else
    display('you are overfitting')
    FIT3D=Data3D;
    FIT3D=zeros(size(Data3D(:,:,:)));
    Residuals=FIT3D;
    Residuals=(Data3D-FIT3D).*Mask;
end;

function model = Create_model(x1,y1,z1,dim,order);
Nsize=[1 4 10 20 35];
N=Nsize(order+1);
model=double(zeros(length(x1),N));
    if order>=0
    model(1:length(x1),1)=1;
        %         model(1:length(Indices),1)=1;
    end
    %first order
    if order>=1
        model(:,2)=reshape(x1-dim(1)/2,length(x1),1);%x -siemens
        model(:,3)=reshape(y1-dim(2)/2,length(x1),1);%y -siemens
        model(:,4)=reshape(z1-dim(3)/2,length(x1),1);%z -siemens
    end
    %second order
    if order>=2
        model(:,5) =model(:,2).*model(:,2);%x^2y^0z^0
        model(:,6) =model(:,2).*model(:,3);%x^1y^1z^0 -siemens
        model(:,7) =model(:,3).*model(:,3);%x^0y^2z^0
        model(:,8) =model(:,3).*model(:,4);%x^0y^1z^1 -siemens
        model(:,9) =model(:,4).*model(:,4);%x^0y^0z^2 -siemens
        model(:,10)=model(:,4).*model(:,2);%x^1y^0z^1 -siemens
    end
    %third order
    if order>=3
        model(:,11) =model(:,5).*model(:,2); %x^3y^0z^0
        model(:,14) =model(:,7).*model(:,3); %x^0y^3z^0
        model(:,18) =model(:,9).*model(:,4); %x^0y^0z^3
        model(:,15) =model(:,8).*model(:,2); %x^1y^1z^1
        model(:,12) =model(:,6).*model(:,2); %x^2y^1z^0
        model(:,13) =model(:,7).*model(:,2); %x^1y^2z^0
        model(:,16) =model(:,9).*model(:,2); %x^1y^0z^2
        model(:,17) =model(:,10).*model(:,2);%x^2y^0z^1
        model(:,19) =model(:,9).*model(:,3); %x^0y^1z^2
        model(:,20) =model(:,8).*model(:,3); %x^0y^2z^1
    end
    if order>=4
        model(:,21) =model(:,11).*model(:,2); %x^4y^0z^0
        model(:,22) =model(:,14).*model(:,3); %x^0y^4z^0
        model(:,23) =model(:,18).*model(:,4); %x^0y^0z^4
        model(:,24) =model(:,11).*model(:,3); %x^3y^1z^0
        model(:,25) =model(:,11).*model(:,4); %x^3y^0z^1
        model(:,26) =model(:,5).*model(:,7); %x^2y^2z^0
        model(:,27) =model(:,5).*model(:,9); %x^2y^0z^2
        model(:,28) =model(:,17).*model(:,3);%x^2y^1z^1
        model(:,29) =model(:,14).*model(:,4); %x^0y^3z^1
        model(:,30) =model(:,14).*model(:,2); %x^1y^3z^0
        model(:,31) =model(:,13).*model(:,4); %x^1y^2z^1
        model(:,32) =model(:,20).*model(:,4); %x^0y^2z^2
        model(:,33) =model(:,18).*model(:,3); %x^0y^1z^3
        model(:,34) =model(:,18).*model(:,2); %x^1y^0z^3
        model(:,35) =model(:,15).*model(:,4); %x^1y^1z^2
    end

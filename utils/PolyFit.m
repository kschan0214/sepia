function [FIT3D,Residuals,b]=PolyFit(Data3D,Mask,order);

dim=size(Data3D);


%% for volume fitting


    Indices = cast(find(Mask(:,:,:)),'like',Data3D);
    [x1,y1,z1]=ind2sub(size(Data3D(:,:,:)),Indices);
    x1 = cast(x1,'like',Data3D);
    y1 = cast(y1,'like',Data3D);
    z1 = cast(z1,'like',Data3D);
    R=Data3D(Indices);
    %%%%%%%%%%this kind of works... but it would be better to make it 3D
    model(1:length(Indices),1) = zeros(length(Indices),1, 'like', Data3D);
%zeroth order
if length(Indices)>(3*order)^2
    if order>=0
        model = ones(length(Indices),1, 'like', Data3D);
%         model(1:length(Indices),1)=1;
    end
        %first order
    if order>=1
        model(:,2)=reshape(x1-dim(1)/2,length(Indices),1);%x -siemens
        model(:,3)=reshape(y1-dim(2)/2,length(Indices),1);%y -siemens
        model(:,4)=reshape(z1-dim(3)/2,length(Indices),1);%z -siemens
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
%         b=model\R;
        b=pinv(model)*R;
        Fit=model*b;
        media=mean(squeeze(Fit));
        desvio=std(squeeze(Fit));
        clear model
        FIT3D=Data3D;
        FIT3D=zeros(size(Data3D(:,:,:)), 'like', Data3D);
        for pos=1:length(x1)
            FIT3D(x1(pos),y1(pos),z1(pos))=Fit(pos);
        end;
        clear Fit
        FIT3D(isnan(FIT3D))=0;
        Residuals=FIT3D;
        Residuals=(Data3D-FIT3D).*Mask;
%     figure
%     subplot(331)
%     imagesc(Data3D(:,:,round(dim(3)/2)),[(media-desvio) (media+desvio)]);colorbar
%     subplot(332)
%     imagesc(FIT3D(:,:,round(dim(3)/2)),[(media-desvio) (media+desvio)]);colorbar
%     subplot(333)
%     imagesc(Residuals(:,:,round(dim(3)/2)));colorbar
%     subplot(337)
%     imagesc(squeeze(Data3D(round(dim(1)/2),:,:)),[(media-desvio) (media+desvio)]);colorbar
%     subplot(338)
%     imagesc(squeeze(FIT3D(round(dim(1)/2),:,:)),[(media-desvio) (media+desvio)]);colorbar
%     subplot(339)
%     imagesc(squeeze(Residuals(round(dim(1)/2),:,:)));colorbar
%     subplot(334)
%     imagesc(squeeze(Data3D(:,round(dim(2)/2),:)),[(media-desvio) (media+desvio)]);colorbar
%     subplot(335)
%     imagesc(squeeze(FIT3D(:,round(dim(2)/2),:)),[(media-desvio) (media+desvio)]);colorbar
%     subplot(336)
%     imagesc(squeeze(Residuals(:,round(dim(2)/2),:)));colorbar
else
        FIT3D=Data3D;
        FIT3D=zeros(size(Data3D(:,:,:)));
        Residuals=FIT3D;
        Residuals=(Data3D-FIT3D).*Mask;
end;
        
% %     figure
% %     keyboard
% %     
% %     imagesc(FIT3D());colorbar

function [fieldMap T2map M0map Confidence]=T2starAndFieldCalc(Magnitude,Phase,te,varargin);
% usage
%[fieldMap T2correctmap M0map Confidence] =T2starAndFieldCalc(Magnitude,Phase,te,varargin);
% if a varagin is not introduced the method used for T2* calculations is
% a simple integral method.
% varargin{1} is the method used to perform the calculation:
% 'fast'- uses the integral of the curve to esimtate T2*
% 'pinv_reweight'- uses the log of the data to perform the line fit, but uses a reweighting of the points to make it sensitive to the lack of confidence on points with smaller intensities or bigger phase rolls
% 'lsqr_reweight'- uses the log of the data to perform the line fit, but uses a reweighting of the points to make it sensitive to the lack of confidence on points with smaller intensities or bigger phase rolls
% 'else'- it simply does the Mathworks recommended exponential fit
%  varargin{2} ==1 a correction of the field gradient is applied exponential decay
% varargin{3} mask to be applied to the data

%%
% output of the fieldMap is in rad/s
% keyboard
method='fast';
dims=size(Magnitude)

if nargin==3
    method='fast';
    gradcorrection=0;
    mask=ones(dims(1:3));
else
    if nargin>=4
        
        if ~isempty(varargin{1})
            method=varargin{1};
        else
            method='fast';
        end;
        if nargin>=5
            
            if ~isempty(varargin{2})
                gradcorrection=varargin{2};
            else
                gradcorrection=0;
            end;
            if nargin>=6
                
                if ~isempty(varargin{3})
                    mask=varargin{3};
                else
                    mask=ones(dims(1:3));
                end;
                
            else
                mask=ones(dims(1:3));
                
            end
            
        else
            gradcorrection=0;
            mask=ones(dims(1:3));
            
        end
        
    else
        method='fast';
        
        gradcorrection=0;
        mask=ones(dims(1:3));
        
    end
end;
%  keyboard;
pos=round(centerofmass(Magnitude))
% pos=round(dims/2); %almost on resonance position



temp3=angle(exp(1i*Phase(:,:,:,2:end))./exp(1i*Phase(:,:,:,1:end-1)));

% unwrapping consecutive echoes

for k=1:dims(4)-1
    k
    
    temp=make_nii(single(((temp3(:,:,:,k)))));
    %     temp.img=temp.img;
    % puts noise outside the mask to give freedom to the software in regions outside the mask;
    temp.img=temp.img.*mask+(1-mask).*(rand(dims(1:3))*2*pi-pi);
    temp.hdr.dime.datatype=16;
    temp.hdr.dime.bitpix=16;
    save_nii(temp,'temp.hdr');
%     unix('sh  /home/rebelo/Documents/MATLAB/phase_unwrapping/test temp.hdr')
    unix('sh  /home/mrphys/kwocha/Tools/phase_unwrap/unwrap temp.hdr')
    temp=load_nii('uwtemp.hdr');
    temp2(1:dims(1),1:dims(2),1:dims(3),k)=temp.img-round(temp.img(pos(1),pos(2),pos(3))/(2*pi))*2*pi;
end;
%
temp3=cumsum(temp2,4);
for k=1:dims(4)-1
    temp3(:,:,:,k)=temp3(:,:,:,k)/(te(k+1)-te(1));
end;
denominator=zeros(dims(1:3));
numerator=zeros(dims(1:3));
for k=1:dims(4)-1
    numerator= numerator + (temp3(:,:,:,k)) ...
        .*abs((te(k+1)-te(1))*Magnitude(:,:,:,k+1)).^2 ...
        ./(abs(Magnitude(:,:,:,k+1)).^2+abs(Magnitude(:,:,:,1)).^2);
    denominator= denominator + ...
        abs((te(k+1)-te(1))*Magnitude(:,:,:,k+1)).^2 ...
        ./(abs(Magnitude(:,:,:,k+1)).^2+abs(Magnitude(:,:,:,1)).^2);
    
end;

fieldMap=numerator./denominator;
fieldMap(isnan(fieldMap))=0;
fieldMap(isinf(fieldMap))=0;
fieldMap=fieldMap.*mask;
%  keyboard

clear numerator;clear denominator;clear temp;clear temp2;clear temp3;
if gradcorrection==1
    display(['calculating field gradients' ])
    
    % corrects the exponential decay using the gradient estimate
    [fieldgradient(:,:,:,1), fieldgradient(:,:,:,2), fieldgradient(:,:,:,3)]=gradient_3D(fieldMap,[],0);
    noisel=norm(fieldgradient(:),1)/numel(fieldgradient)
    fieldMapdn=Wavedec3Denoising(fieldMap,noisel,4,'sym4','verysoft');
    fieldMapdn=real(fieldMapdn);
    % image_view4([gradient_3D(fieldMapdn),gradient_3D(fieldMap)]);
    
    [fieldgradient(:,:,:,1), fieldgradient(:,:,:,2), fieldgradient(:,:,:,3)]=gradient_3D(fieldMapdn,[],0);
    slicefraction=1;
    for t=1:length(te)
        weight(:,:,:,t)=abs(sinc(fieldgradient(:,:,:,1)*te(t)/2/pi)).*...
            abs(sinc(fieldgradient(:,:,:,2)*te(t)/2/pi)).*...
            abs(sinc(fieldgradient(:,:,:,3)*te(t)/2/pi));
    end;
    Magnitude=Magnitude./weight;
else
    weight=ones(dims);
    fieldgradient=zeros(dims(1:3));
end;
Magnitude(isinf(Magnitude))=0;
Magnitude(isnan(Magnitude))=0;
display(['calculate T2* map' ])

if strcmp(method,'fast')
    temp=0;
    for k=1:dims(4)-1
        temp=temp+0.5*(Magnitude(:,:,:,k)+Magnitude(:,:,:,k+1))*(te(k+1)-te(k));
    end;
    % very fast estimation
    T2map=temp./(Magnitude(:,:,:,1)-Magnitude(:,:,:,end));
    T2map(isinf(T2map))=min(te)/20;
    T2map((T2map)<0)=max(te)*20;
    T2map((T2map)>max(te)*20)=max(te)*20;
    T2map(isnan(T2map))=max(te)*20;
    M0map=(Magnitude(:,:,:,1).*exp(te(1)./T2map));M0map(isinf(M0map))=0;M0map(isnan(M0map))=0;
    %     M0=reshape(Magnitude,[prod(dims(1:3)) dims(4)])*pinv(exp(-1./T2map(:)*te));
    % elseif and(gradcorrection==1,strcmp(method,'pinv_reweighted'))
elseif (strcmp(method,'pinv_reweighted'))
    % it is a very slow method
    temp=0;
    for k=1:dims(4)-1
        temp=temp+0.5*(Magnitude(:,:,:,k)+Magnitude(:,:,:,k+1))*(te(k+1)-te(k));
    end;
    T2map=temp./(Magnitude(:,:,:,1)-Magnitude(:,:,:,end));
        T2map(isinf(T2map))=min(te)/20;
    T2map((T2map)<0)=max(te)*20;
    T2map((T2map)>max(te)*20)=max(te)*20;
    T2map(isnan(T2map))=max(te)*20;

    M0map=(Magnitude(:,:,:,1).*exp(te(1)./T2map));M0map(isinf(M0map))=0;M0map(isnan(M0map))=0;
    % very fast estimation done
    'fast estimation done'
    mask(sum(abs(fieldgradient)*te(1)/(2*pi),4)>0.7)=0;%masks data when the fieldmap gradient means that by te(1) there was already one full turn of the spins
    mask(sum(abs(fieldgradient)*te(end)/(2*pi),4)<0.4)=0;%masks data when the fieldmap gradient only creates a decay of 25% by the last te
    %% I AM HERE MAKING IT ONLY FOR THE pixels that need to be corrected
    indices=find(mask==1);
    [x , y ,z]=ind2sub(dims,indices);
    tic
    te2=zeros(2,length(te));gre_5Echo0p8mmres
    for k=1:length(indices)
        %     k
        te2(1,:)=weight(x(k),y(k),z(k),:);
        te2(2,:)=reshape(te,[1,length(te)]).*reshape(weight(x(k),y(k),z(k),:),[1,length(te)]);
        BallLog=squeeze(log(Magnitude(x(k),y(k),z(k),:)).*weight(x(k),y(k),z(k),:))';
        R2line=-(BallLog*pinv(te2));
        
        T2map(x(k),y(k),z(k))=1./R2line(2);
        M0map(x(k),y(k),z(k))=exp(-R2line(1));
    end;
    toc
    
    % elseif and(gradcorrection==1,strcmp(method,'lsqr_reweight'))
elseif (strcmp(method,'lsqr_reweight'))
    % starts by doing the fast method
    temp=0;
    for k=1:dims(4)-1
        temp=temp+0.5*(Magnitude(:,:,:,k)+Magnitude(:,:,:,k+1))*(te(k+1)-te(k));
    end;gre_5Echo0p8mmres
    T2map=temp./(Magnitude(:,:,:,1)-Magnitude(:,:,:,end));
    T2map(isinf(T2map))=min(te)/20;
    T2map((T2map)<0)=max(te)*20;
    T2map((T2map)>max(te)*20)=max(te)*20;
    T2map(isnan(T2map))=max(te)*20;

    M0map=log(Magnitude(:,:,:,1).*exp(te(1)./T2map));M0map(isinf(M0map))=0;M0map(isnan(M0map))=0;
    M0maptemp=M0map;
    'fast estimation done'
    % very fast estimation done
    
    weights=Magnitude.*weight;
    % log of the decay
    Y=log(Magnitude);
    
    %removing useless points that jsut create problems
    weights(isinf(Y))=0;mask(prod(weights,4)==0)=0;
    Y(isinf(Y))=0;mask(prod(Y,4)==0)=0;
    X0=cat(4,M0maptemp,abs(1./T2map));
    % X0(X0<0)=0;
    X0=X0(:);
    X0(isnan(X0))=0;
    
%   keyboard;
   
    %only does the fitting if
    mask(sum(abs(fieldgradient)*te(1)/(2*pi),4)>0.7)=0;%masks data when the fieldmap gradient means that by te(1) there was already one full turn of the spins
    mask(sum(abs(fieldgradient)*te(end)/(2*pi),4)<0.3)=0;%masks data when the fieldmap gradient only creates a decay of 25% by the last te
    weights=bsxfun(@times,weights,mask);
    
%     keyboard
    X = lsqr(@(x,tflag)T2starfitlsqr(x,te,weights(:),tflag),Y(:).*weights(:),1e-3,100,[],[],X0);
    
    %     norm(T2starfitlsqr(X0,te,weights(:),'notransp')-Y(:).*weights(:))
    %     norm(T2starfitlsqr(X,te,weights(:),'notransp')-Y(:).*weights(:))
    
    X=reshape(X,[size(Magnitude,1),size(Magnitude,2),size(Magnitude,3),2]);
    %     X0=reshape(X0,[size(Magnitude,1),size(Magnitude,2),size(Magnitude,3),2]);
    % image_view3(sum(reshape((T2starfitlsqr(X,te,weights(:),'notransp')-Y(:).*weights(:)).^2,[size(Magnitude)]),4))
    M0map= exp(X(:,:,:,1)); M0map(mask==0)=exp(M0maptemp(mask==0));
    R2map= X(:,:,:,2); R2map(mask==0)=1./T2map(mask==0);
    T2map=1./R2map;
    
else
    %doing an exponential fit
    if (prod(dims)*16<75e6)
        volmesliceline=1;
    elseif (prod(dims([1,2,4]))*24<50e6)
        volmesliceline=2;
    else
        volmesliceline=3;
    end;
    [M0map,T2map]=T2mapNifti_Mathworld(Magnitude,te,volmesliceline);
end;

T2map(isinf(T2map))=min(te)/20;
T2map((T2map)<0)=max(te)*20;
T2map((T2map)>max(te)*20)=max(te)*20;
T2map(isnan(T2map))=max(te)*20;

reconData=weight.*bsxfun(@times,M0map,exp(bsxfun(@times,(-1./T2map+i*fieldMap),reshape(te,[1 1 1 dims(4)]))));
OriData=Magnitude.*weight.*exp(i*Phase);
reconData=(bsxfun(@times,reconData,exp(-i*angle(reconData(:,:,:,1)))));
OriData=(bsxfun(@times,OriData,exp(-i*angle(OriData(:,:,:,1)))));
% Confidence=sum((bsxfun(@times,reconData,exp(-i*angle(reconData)))-bsxfun(@times,OriData,exp(-i*angle(OriData)))).^2,4);%./sum(abs(OriData).^2,4);
Confidence=sum(abs(reconData-OriData).^2,4)./sum(abs(OriData).^2,4);
% Confidenceb=sum((bsxfun(@times,reconData,exp(-i*angle(reconData)))-bsxfun(@times,OriData,exp(-i*angle(OriData)))).^2,4)./sum(abs(OriData).^2,4);


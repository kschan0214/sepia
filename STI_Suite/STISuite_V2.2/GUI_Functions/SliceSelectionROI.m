function SliceSelectionROI(hObject,handles)
ImageAxis=handles.ImageAxisNum;
ss=size(handles.Mi,ImageAxis);

set(handles.FrameSelect, 'Max', ss,'SliderStep',[1/ss 2/ss]); 
if get(handles.FrameSelect, 'value')>ss;
    set(handles.FrameSelect, 'value', ss/2); 
end

value = get(handles.FrameSelect,'Value');
if value<1;
    value=1;
end
DataSetNumber=round(str2double(get(handles.CoilNum,'String')));

if get(handles.Scaling,'Value')==1
    scalling='on';
else
    scalling='off';
end


%% 3D display
switch ImageAxis
    case 1
        Image2=squeeze(handles.Mi(value,:,:,DataSetNumber));
        AlphaTemp=squeeze(handles.M2(value,:,:));
    case 2
        Image2=squeeze(handles.Mi(:,value,:,DataSetNumber));
        AlphaTemp=squeeze(handles.M2(:,value,:));
    case 3
        Image2=squeeze(handles.Mi(:,:,value,DataSetNumber));
        AlphaTemp=squeeze(handles.M2(:,:,value));
end


AlphaTempValue=round(AlphaTemp*100)-handles.baseAlpha*100;

img1r=(AlphaTempValue==1|AlphaTempValue==4|AlphaTempValue==5);
img1g=(AlphaTempValue==2|AlphaTempValue==4|AlphaTempValue==6);
img1b=(AlphaTempValue==3|AlphaTempValue==5|AlphaTempValue==6);

img1r=img1r+(AlphaTempValue==7)*0.5;
img1g=img1g+(AlphaTempValue==8)*0.5;
img1b=img1b+(AlphaTempValue==9)*0.5;

img1=zeros([size(img1r),3]);

img1(:,:,1)=img1r;
img1(:,:,2)=img1g;
img1(:,:,3)=img1b;

if get(handles.AbsButton,'value')==1
    Image2=abs(Image2);
else
    Image2=real(Image2);
end

imgrgb=img1r+img1g+img1b;

if get(handles.FlipUpDown,'value')==1
    Image2=imrotate(Image2,180);
    img1=imrotate(img1,180);
    imgrgb=imrotate(imgrgb,180);
    
end
if get(handles.Rotate90,'value')==1
    Image2=imrotate(Image2,90);
    img1=imrotate(img1,90);
    imgrgb=imrotate(imgrgb,90);
end

Image2=permute(Image2(:,:),[2 1]);
img1=permute(img1(:,:,:),[2 1 3]);
imgrgb=permute(imgrgb(:,:),[2 1]);

alphavalue=ones(size(imgrgb));
set(handles.himage2,'CData',Image2,'parent',handles.ImageAxes2);
set(handles.himage1,'Cdata',img1);

alphavalue(imgrgb>0)=get(handles.AlphaValueSlider, 'value');
set(handles.himage2,'AlphaData',squeeze(alphavalue));

if strcmp(scalling, 'on')
    set(handles.ImageAxes,'XLim',[0 size(img1,2)],'YLim',[0 size(img1,1)]);
    set(handles.ImageAxes2,'XLim',[0 size(img1,2)],'YLim',[0 size(img1,1)]);
end

STR=get(handles.LeftResolution,'string');
eval(['resolution=' STR ';']);

res0=(resolution/resolution(2));

switch handles.ImageAxisNum
case 1
    res=[ res0(3) res0(2) res0(1) ];
    if get(handles.Rotate90,'value')==1
        res=[res0(3) res0(2) res0(1) ];
    end
    set (handles.ImageAxes2,'DataAspectRatio',res);
    set (handles.ImageAxes,'DataAspectRatio',res);
case 2
    res=[res0(3) res0(1) res0(2) ];
    if get(handles.Rotate90,'value')==1
        res=[res0(3) res0(1) res0(2) ];
    end
    set (handles.ImageAxes2,'DataAspectRatio',res);
    set (handles.ImageAxes,'DataAspectRatio',res);
case 3
    res=[ res0(2) res0(1) res0(3) ];
    if get(handles.Rotate90,'value')==1
        res=[res0(2) res0(1) res0(3) ];
    end
    set (handles.ImageAxes2,'DataAspectRatio',res);
    set (handles.ImageAxes,'DataAspectRatio',res);
end

end




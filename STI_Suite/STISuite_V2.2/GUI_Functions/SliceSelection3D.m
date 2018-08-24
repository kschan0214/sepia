function SliceSelection3D(hObject,handles)
ImageAxis=handles.ImageAxisNum;
ss=size(handles.Mi,ImageAxis);
set(handles.FrameSelect3D, 'Max', ss,'SliderStep',[1/ss 2/ss]); 


if get(handles.FrameSelect3D, 'value')>ss;
    set(handles.FrameSelect3D, 'value', ss/2); 
end

value = round(get(handles.FrameSelect3D,'Value'));
if value<1;
    value=1;
end

if get(handles.AbsButton,'value')==1
    switch ImageAxis
        case 1
            Image=real(squeeze(handles.Mi(value,:,:,:)));
        case 2
            Image=real(squeeze(handles.Mi(:,value,:,:)));
        case 3
             Image=real(squeeze(handles.Mi(:,:,value,:)));
    end
else
     switch ImageAxis
        case 1
            
            Image=real(squeeze(handles.Mi(value,:,:,:)));
        case 2
            Image=real(squeeze(handles.Mi(:,value,:,:)));
        case 3
            Image=real(squeeze(handles.Mi(:,:,value,:)));
    end
end


handles.LeftImax1=str2double(get(handles.RGBhigh,'string'));
handles.LeftImin1=str2double(get(handles.RGBlow,'string'));

Image=(Image-handles.LeftImin)/(handles.LeftImax-handles.LeftImin);
Image=(Image-handles.LeftImin1)/(handles.LeftImax1-handles.LeftImin1);
Image(Image>1)=1;
Image(Image<0)=0;


if get(handles.FlipUpDown,'value')==1
    Image=imrotate(Image,180);
end
if get(handles.Rotate90,'value')==1
    Image=imrotate(Image,90);
end

Image=permute(Image(:,end:-1:1,:),[2 1 3]);

set(handles.himage1,'CData',squeeze(Image),'parent',handles.ImageAxes);
set(handles.ImageAxes,'XLim',[0 size(Image,2)],'YLim',[0 size(Image,1)]);


%%
rightimage=0;
try

    if get(handles.AbsButton,'value')==1
    switch ImageAxis
        case 1
            Image=real(squeeze(handles.M2(value,:,:,:)));
        case 2
            Image=real(squeeze(handles.M2(:,value,:,:)));
        case 3
            Image=real(squeeze(handles.M2(:,:,value,:)));
    end
else
     switch ImageAxis
        case 1
            Image=real(squeeze(handles.M2(value,:,:,:)));
        case 2
            Image=real(squeeze(handles.M2(:,value,:,:)));
        case 3
            Image=real(squeeze(handles.M2(:,:,value,:)));
    end
end


handles.RightImax1=str2double(get(handles.RightHigh,'string'));
handles.RightImin1=str2double(get(handles.RightLow,'string'));

Image=(Image-handles.RightImin)/(handles.RightImax-handles.RightImin);
Image=(Image-handles.RightImin1)/(handles.RightImax1-handles.RightImin1);
Image(Image>1)=1;
Image(Image<0)=0;


if get(handles.FlipUpDown,'value')==1
    Image=imrotate(Image,180);
end
if get(handles.Rotate90,'value')==1
    Image=imrotate(Image,90);
end

Image=permute(Image(:,end:-1:1,:),[2 1 3]);
set(handles.himage2,'CData',squeeze(Image),'parent',handles.ImageAxes2);
set(handles.ImageAxes2,'XLim',[0 size(Image,2)],'YLim',[0 size(Image,1)]);

rightimage=1;
end



STR=get(handles.LeftResolution,'string');
eval(['resolution=' STR ';']);

res0=(resolution/resolution(2));
switch ImageAxis
case 1
    res=[ res0(3) res0(2) res0(1) ];
    if get(handles.Rotate90,'value')==1
        res=[res0(3) res0(2) res0(1) ];
    end
    set (handles.ImageAxes,'DataAspectRatio',res);
    if rightimage
    set (handles.ImageAxes2,'DataAspectRatio',res);
    end
case 2
    res=[ res0(3) res0(1) res0(2) ];
    if get(handles.Rotate90,'value')==1
        res=[res0(3) res0(1) res0(2) ];
    end
    set (handles.ImageAxes,'DataAspectRatio',res);
    if rightimage
    set (handles.ImageAxes2,'DataAspectRatio',res);
    end
case 3
    res=[ res0(2) res0(1) res0(3) ];
    if get(handles.Rotate90,'value')==1
        res=[res0(2) res0(1) res0(3) ];
    end
    set (handles.ImageAxes,'DataAspectRatio',res);
    if rightimage
    set (handles.ImageAxes2,'DataAspectRatio',res);
    end
end


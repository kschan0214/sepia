function SliceSelection(hObject,handles)

ImageAxis=handles.ImageAxisNum;
ss=size(handles.Mi,ImageAxis);
set(handles.FrameSelect, 'Max', ss,'SliderStep',[1/ss 2/ss]); 
if get(handles.FrameSelect, 'value')>ss;
    set(handles.FrameSelect, 'value', ss/2); 
end

value = floor( get(handles.FrameSelect,'Value'));
if value<1;
    value=1;
end

try
if get(handles.Scaling,'Value')==1
    scalling='on';
else
    scalling='off';
end


catch
    scalling='on';
end



value4 = get(handles.Coilslider,'Value');

if value4 ==0
    return
end

if get(handles.AbsButton,'value')==1
    switch ImageAxis
        case 1
            Image=abs(squeeze(handles.Mi(value,:,:,value4)));
        case 2
            Image=abs(squeeze(handles.Mi(:,value,:,value4)));
        case 3
             Image=abs(squeeze(handles.Mi(:,:,value,value4)));
    end
else
     switch ImageAxis
        case 1
            Image=real(squeeze(handles.Mi(value,:,:,value4)));
        case 2
            Image=real(squeeze(handles.Mi(:,value,:,value4)));
        case 3
            Image=real(squeeze(handles.Mi(:,:,value,value4)));
    end
end


if get(handles.FlipUpDown,'value')==1
    Image=Image(end:-1:1,:);
end
if get(handles.Rotate90,'value')==1
    Image=permute(Image,[2 1]);
end

Image=permute(Image(:,end:-1:1),[2 1]);
set(handles.himage1,'CData',Image,'parent',handles.ImageAxes);
if strcmp(scalling, 'on')
set(handles.ImageAxes,'XLim',[0 size(Image,2)],'YLim',[0 size(Image,1)]);
end
rightimage=0;


try
if get(handles.AbsButton,'value')==1
    switch ImageAxis
        case 1
            Image2=abs(squeeze(handles.M2(value,:,:,value4)));
        case 2
            Image2=abs(squeeze(handles.M2(:,value,:,value4)));
        case 3
             Image2=abs(squeeze(handles.M2(:,:,value,value4)));
    end
else
     switch ImageAxis
        case 1
            Image2=real(squeeze(handles.M2(value,:,:,value4)));
        case 2
            Image2=real(squeeze(handles.M2(:,value,:,value4)));
        case 3
            Image2=real(squeeze(handles.M2(:,:,value,value4)));
    end
end
if get(handles.FlipUpDown,'value')==1
    Image2=Image2(end:-1:1,:);
end
if get(handles.Rotate90,'value')==1
    Image2=permute(Image2,[2 1]);
end

Image2=permute(Image2(:,end:-1:1),[2 1]);
set(handles.himage2,'CData',Image2,'parent',handles.ImageAxes2);
if strcmp(scalling, 'on')
set(handles.ImageAxes2,'XLim',[0 size(Image,2)],'YLim',[0 size(Image,1)]);
end
rightimage=1;
catch
    try
    if get(handles.AbsButton,'value')==1
        switch ImageAxis
            case 1
                Image2=abs(squeeze(handles.M2(value,:,:,1)));
            case 2
                Image2=abs(squeeze(handles.M2(:,value,:,1)));
            case 3
                 Image2=abs(squeeze(handles.M2(:,:,value,1)));
        end
    else
         switch ImageAxis
            case 1
                Image2=real(squeeze(handles.M2(value,:,:,1)));
            case 2
                Image2=real(squeeze(handles.M2(:,value,:,1)));
            case 3
                Image2=real(squeeze(handles.M2(:,:,value,1)));
        end
    end
    if get(handles.FlipUpDown,'value')==1
        Image2=Image2(end:-1:1,:);
    end
    if get(handles.Rotate90,'value')==1
        Image2=permute(Image2,[2 1]);
    end

    Image2=permute(Image2(:,end:-1:1),[2 1]);
    set(handles.himage2,'CData',Image2,'parent',handles.ImageAxes2);
    if strcmp(scalling, 'on')
    set(handles.ImageAxes2,'XLim',[0 size(Image,2)],'YLim',[0 size(Image,1)]);
    end
    rightimage=1;    
    end
end
%%


STR=get(handles.LeftResolution,'string');
eval(['resolution=' STR ';']);
res0=(resolution/resolution(1));
switch ImageAxis
case 1
    res=[res0(3) res0(2)  res0(1) ];
    if get(handles.Rotate90,'value')==1
        res=[res0(2) res0(3) res0(1) ];
    end
    set (handles.ImageAxes,'DataAspectRatio',res);
    if rightimage
        set (handles.ImageAxes2,'DataAspectRatio',res);
    end
case 2
    res=[ res0(3) res0(1) res0(2) ];
    if get(handles.Rotate90,'value')==1
        res=[res0(1) res0(3) res0(1) ];
    end
    set (handles.ImageAxes,'DataAspectRatio',res);
    if rightimage
        set (handles.ImageAxes2,'DataAspectRatio',res);
    end
case 3
    res=[res0(2) res0(1)  res0(3) ];
    if get(handles.Rotate90,'value')==1
        res=[res0(1) res0(2) res0(3) ];
    end
    set (handles.ImageAxes,'DataAspectRatio',res);
    if rightimage
        set (handles.ImageAxes2,'DataAspectRatio',res);
    end
end



function VisualizeROI(hObject,handles,ctrl)
try
    switch ctrl
        case 1
            myvarname=get(handles.VariableSelect,'string');
            handles.Mi=evalin('base', myvarname);
            ss=size(handles.Mi);
            SS3=1;
        case 2
            myvarname=get(handles.AddMoreImages,'string');
            MoreImage=evalin('base', myvarname);
            clc
            SS1=size(handles.Mi,4);
            SS2=size(MoreImage,4);
            handles.Mi(:,:,:,SS1+(1:SS2))=MoreImage;
            SS3=SS1+1;
            ss=size(handles.Mi);
    end
catch 
    set(handles.status,'string','Error! No Image Selected ...')
    return
end
if length(ss)<3; ss(3)=1; end
set(handles.FrameSelect, 'Max', ss(3),'SliderStep',[1/ss(3) 2/ss(3)]); 
set(handles.FrameSelect, 'visible','on'); 

handles.M2=ones(size(handles.Mi,1),size(handles.Mi,2),size(handles.Mi,3));


handles.img1(:,:,1)=real(handles.M2(:,:,round(ss(3)/2)));
handles.img1(:,:,3)=0;


Image=handles.img1;
Image=permute(Image,[2 1 3]);
handles.himage1=imshow(Image,'parent',handles.ImageAxes);

Image2=real(handles.Mi(:,:,round(ss(3)/2),SS3));
Image2=permute(Image2(:,:),[2 1]);
handles.himage2=imshow(Image2,[],'parent',handles.ImageAxes2);


AlphaMap=handles.M2(:,:,round(ss(3)/2));
AlphaMap=(AlphaMap>0);

AlphaMap=permute(AlphaMap(:,end:-1:1),[2 1]);

set(handles.himage2,'AlphaData',AlphaMap);
set(handles.himage1,'CDataMapping','scaled')
if get(handles.Grey,'value')
    colormap gray
else 
    colormap jet
end

sss=size(handles.Mi,4);

if sss>1
    set(handles.Coilslider,'Max', sss,'SliderStep',[1/(sss-1) 2/(sss-1)]); 
    if get(handles.Coilslider, 'value')>sss;
        set(handles.Coilslider, 'value', sss/2); 
    end
    set(handles.Coilslider,'value', 1); 
end

Xlim2=get (handles.ImageAxes2,'Xlim');
Ylim2=get (handles.ImageAxes2,'ylim');

XYLimMax=max([Xlim2(2)-Xlim2(1) Ylim2(2)-Ylim2(1)]);

set (handles.ImageAxes2,'Xlim',[0 XYLimMax]+0.5,'Ylim',[0 XYLimMax]+0.5);
set (handles.ImageAxes,'Xlim',[0 XYLimMax]+0.5,'Ylim',[0 XYLimMax]+0.5);
set(handles.FrameSelect, 'value', round(ss(3)/2)); 
set(handles.status,'string','Images loaded.')
set(handles.CoilNum,'String',num2str(1));
set (handles.ImageAxes2,'ydir','normal');
set (handles.ImageAxes,'ydir','normal');

handles.ImageAxisNum=3;
handles.NumberSelect=1;
handles.baseAlpha=0.5;

guidata(hObject, handles);

set(handles.ImageAxes,'XLim',[0 size(Image,2)],'YLim',[0 size(Image,1)]);
set(handles.ImageAxes2,'XLim',[0 size(Image,2)],'YLim',[0 size(Image,1)]);

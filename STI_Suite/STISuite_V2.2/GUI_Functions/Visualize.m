function Visualize(hObject,handles,ctrl)
try
    switch ctrl
        case 1
            myvarname=get(handles.VariableSelect,'string');
            handles.Mi=evalin('base', myvarname);
            ss=size(handles.Mi);
        case 2
            myvarname=get(handles.RightImageLoad,'string');
            handles.M2=evalin('base', myvarname);
            ss=size(handles.M2);
        case 3
            myvarname=get(handles.ProcessedPhaseName,'string');
            handles.M2=evalin('base', myvarname);
            ss=size(handles.M2);
        case 4
            myvarname=get(handles.NewMaskName,'string');
            handles.M2=evalin('base', myvarname);
            ss=size(handles.M2);
        case 5
            myvarname=get(handles.unwrappedphase,'string');
            handles.M2=evalin('base', myvarname);
            ss=size(handles.M2);
        case 6
            myvarname=get(handles.BrainMaskName,'string');
            handles.M2=evalin('base', myvarname);
            ss=size(handles.M2);
        case 7
            myvarname=get(handles.loadniiname1,'string');
            handles.Mi=evalin('base', myvarname);
            ss=size(handles.Mi);
        case 8
            myvarname=get(handles.loadniiname2,'string');
            handles.M2=evalin('base', myvarname);
            ss=size(handles.M2);
    end
catch 
    set(handles.status,'string','Error! No Image Selected ...')
    return
end

if length(ss)<3; ss(3)=1; end
set(handles.FrameSelect, 'Max', ss(3),'SliderStep',[1/ss(3) 2/ss(3)]); 
set(handles.FrameSelect, 'visible','on'); 


switch ctrl
    case 1
        Image=real(handles.Mi(:,:,round(ss(3)/2)));
        Image=permute(Image(:,end:-1:1),[2 1]);
        handles.himage1=imshow(Image,[],'parent',handles.ImageAxes);
        
    case 2
        Image=real(handles.M2(:,:,round(ss(3)/2)));
        Image=permute(Image(:,end:-1:1),[2 1]);
        handles.himage2=imshow(Image,[],'parent',handles.ImageAxes2);
    case 3
        Image=real(handles.M2(:,:,round(ss(3)/2)));
        Image=permute(Image(:,end:-1:1),[2 1]);
        handles.himage2=imshow(Image,[],'parent',handles.ImageAxes2);
    case 4
        Image=real(handles.M2(:,:,round(ss(3)/2)));
        Image=permute(Image(:,end:-1:1),[2 1]);
        handles.himage2=imshow(Image,[],'parent',handles.ImageAxes2);
    case 5
        Image=real(handles.M2(:,:,round(ss(3)/2)));
        Image=permute(Image(:,end:-1:1),[2 1]);
        handles.himage2=imshow(Image,[],'parent',handles.ImageAxes2);
    case 6
        Image=real(handles.M2(:,:,round(ss(3)/2)));
        Image=permute(Image(:,end:-1:1),[2 1]);
        handles.himage2=imshow(Image,[],'parent',handles.ImageAxes2);
    case 7
        Image=real(handles.Mi(:,:,round(ss(3)/2)));
        Image=permute(Image(:,end:-1:1),[2 1]);
        handles.himage1=imshow(Image,[],'parent',handles.ImageAxes);
    case 8
        Image=real(handles.M2(:,:,round(ss(3)/2)));
        Image=permute(Image(:,end:-1:1),[2 1]);
        handles.himage2=imshow(Image,[],'parent',handles.ImageAxes2);
end
set(handles.himage1,'CDataMapping','scaled')

if get(handles.Grey,'value')
    colormap gray
else 
    colormap jet
end

sss=size(handles.Mi,4);
if sss>1
    set(handles.Coilslider,'Max', sss,'SliderStep',[1/sss 2/sss]); 
    if get(handles.Coilslider, 'value')>sss;
        set(handles.Coilslider, 'value', sss/2); 
        
    end
    set(handles.Coilslider,'value', 1); 
end
handles.ImageAxisNum=3;
set(handles.FrameSelect, 'value', round(ss(3)/2)); 
set(handles.status,'string','Images loaded.')
guidata(hObject, handles);

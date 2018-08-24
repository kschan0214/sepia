function Brain_ROIs(hObject,handles)

set(handles.himage2,'ButtonDownFcn',@MovePoint1);
set(gcbf,'WindowButtonUpFcn',@stopmove)

function MovePoint1(varargin)
handles = guidata(gcbo);
value = get(handles.FrameSelect,'Value');
switch handles.ImageAxisNum
    case 1
        AlphaTemp=squeeze(handles.M2(value,:,:));
    case 2
        AlphaTemp=squeeze(handles.M2(:,value,:));
    case 3
        AlphaTemp=squeeze(handles.M2(:,:,value));
end
AlphaTemp=permute(AlphaTemp,[2 1]);
handles.AlphaTemp=AlphaTemp;

pt = get(gca,'CurrentPoint');
current=round(pt(1,1:2));
if strcmp(get(gcf,'SelectionType'),'normal')
    handles.AlphaTemp(current(2),current(1))=handles.NumberSelect*0.01+handles.baseAlpha;
else
    handles.AlphaTemp(current(2),current(1))=1;
end
set(handles.himage2,'AlphaData',handles.AlphaTemp);
handles.ColorSelection={[1 0 0],[0 1 0],[0 0 1],[1 1 0],[1 0 1],[0 1 1],[0.5 0 0],[0 0.5 0],[0 0 0.5]};

AlphaTempValue=round(handles.AlphaTemp*100)-handles.baseAlpha*100;
img1r=(AlphaTempValue==1|AlphaTempValue==4|AlphaTempValue==5);
img1g=(AlphaTempValue==2|AlphaTempValue==4|AlphaTempValue==6);
img1b=(AlphaTempValue==3|AlphaTempValue==5|AlphaTempValue==6);

img1r=img1r+(AlphaTempValue==7)*0.5;
img1g=img1g+(AlphaTempValue==8)*0.5;
img1b=img1b+(AlphaTempValue==9)*0.5;


handles.img1=zeros([size(img1r),3]);
handles.img1(:,:,1)=img1r;
handles.img1(:,:,2)=img1g;
handles.img1(:,:,3)=img1b;
img1=handles.img1;
set(handles.himage1,'Cdata',img1);

guidata(gcbo,handles);
set(gcbf,'WindowButtonMotionFcn',@startmove)


function stopmove(varargin)
set(gcf,'Pointer','arrow')
set(gcbf,'WindowButtonMotionFcn','')
handles = guidata(gcbo);
value = get(handles.FrameSelect,'Value');

AlphaTemp=permute(handles.AlphaTemp,[2 1 3]);

if handles.savedata
    switch handles.ImageAxisNum
        case 1
            handles.M2(value,:,:)=AlphaTemp;
        case 2
            handles.M2(:,value,:)=AlphaTemp;
        case 3
            handles.M2(:,:,value)=AlphaTemp;
    end
end
handles.savedata=0;
guidata(gcbo,handles);

function startmove(varargin)
handles = guidata(gcbo);
pt = get(gca,'CurrentPoint');
current=round(pt(1,1:2));
set(gcf,'Pointer','hand');
radius0=str2double(get(handles.radius0,'String'));
SS=size(handles.AlphaTemp);
for xx=(current(2)-radius0):(current(2)+radius0)
    for yy=(current(1)-radius0):(current(1)+radius0)
        if xx>0 && xx <SS(1)+1 && yy>0 && yy <SS(2)+1
          if ((xx-current(2))^2+(yy-current(1))^2)<radius0^2
            if strcmp(get(gcf,'SelectionType'),'normal')
                handles.AlphaTemp(xx,yy)=handles.NumberSelect*0.01+handles.baseAlpha;
                handles.img1(xx,yy,:)=handles.ColorSelection{handles.NumberSelect};
            else
                handles.AlphaTemp(xx,yy)=1;
            end
            
          end
        end
    end
end

set(handles.himage2,'AlphaData',handles.AlphaTemp);
set(handles.himage1,'Cdata',handles.img1);
handles.savedata=1;
guidata(gcbo,handles);







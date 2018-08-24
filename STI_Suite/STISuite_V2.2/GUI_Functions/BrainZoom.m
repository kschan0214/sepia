function BrainZoom(hObject,handles)

Xlim2=get (handles.ImageAxes2,'Xlim');
Ylim2=get (handles.ImageAxes2,'ylim');
set (handles.ImageAxes,'Xlim',Xlim2,'Ylim',Ylim2);
guidata(gcbo,handles);

function ShowMagnitude(hObject,handles,control)
% Wei Li, Duke University, June 2010

switch control
    case 'Mag'
        handles.Mi=handles.magni;
    case 'Chi'
        handles.Mi=handles.susceptibility;
    case 'R2*'
        handles.Mi=handles.R2Star;
end

ss=size(handles.Mi);     
sss=size(handles.Mi,4);
set(handles.himage2,'CData',real(handles.Mi(:,:,round(ss(3)/2))),'parent',handles.ImageAxes2);

if sss>1
set(handles.Coilslider,'Max', sss,'SliderStep',[1/sss 2/sss]); 
if get(handles.Coilslider, 'value')>sss;
    set(handles.Coilslider, 'value', sss/2); 
end
set(handles.Coilslider,'value', 1); 
end
guidata(hObject, handles);
SliceSelectionROI(hObject,handles);

switch control
    case 'Mag'
        A=get(handles.himage2,'CData');
        s1=[min(A(:)) max(A(:))];
    case 'Chi'
        s1=[-1 1]*0.3;
    case 'R2*'
        s1=[0 100];
end

caxis(handles.ImageAxes2,s1);


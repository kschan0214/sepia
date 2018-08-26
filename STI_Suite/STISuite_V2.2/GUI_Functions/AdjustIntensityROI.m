function AdjustIntensityROI(hObject,handles,control)
% Wei Li, Duke University, June 2010

switch control
    case 7
        imcontrast(handles.himage2);
        %get(handles.ImageAxes2,'clim')
    case 3
        A2=get(handles.himage2,'CData');
        s2=[min(A2(:)) max(A2(:))];
        caxis(handles.ImageAxes2,s2);
end



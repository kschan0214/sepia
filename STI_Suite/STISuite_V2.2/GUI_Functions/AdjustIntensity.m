function AdjustIntensity(hObject,handles,control)
% Wei Li, Duke University, June 2010
if control==5
    imcontrast(handles.himage1);
    return
end

if control==6
    A=get(handles.ImageAxes,'Clim');
    set(handles.ImageAxes2,'Clim',A);
    return
end
if control==7
    imcontrast(handles.himage2);
    return
end
if control==8
    A=get(handles.ImageAxes2,'Clim');
    set(handles.ImageAxes,'Clim',A);
    return
end

if control==3
    A=get(handles.himage1,'CData');
    s1=[min(A(:)) max(A(:))];
    try
    A2=get(handles.himage2,'CData');
    s2=[min(A2(:)) max(A2(:))];
    end

end

if control==20
    s1=[-1 1]*0.2;
    s2=[-1 1]*0.2;
else
try
    B=1.02^(get(handles.Intensity,'value')-20);
    C=0.01*(get(handles.Contrast,'value')-20);
    s1=[-0 1]*B+C;
    s2=[-0 1]*B+C;
    end
end

try
caxis(handles.ImageAxes,s1);
caxis(handles.ImageAxes2,s2);
end


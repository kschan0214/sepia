function ColorScheme(hObject,handles)

if get(handles.Grey,'value')
    colormap gray
else 
    colormap jet
end

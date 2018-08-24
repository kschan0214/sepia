function Assign2Base(hObject,handles,ctrl)
switch ctrl
    case 1
        Imag=handles.Mi;
        myvarname=get(handles.Assign2Base,'string');
    case 2
        Imag=handles.M2;
        myvarname=get(handles.RightImageSave,'string');
end

assignin('base',myvarname,double(Imag));

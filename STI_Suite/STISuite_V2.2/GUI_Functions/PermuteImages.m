function PermuteImages(hObject,handles,control)
% Wei Li, 06/14/2014

UniqueID=get(handles.UniqueID,'string');
switch control
    case 1
        %varname=[UniqueID '_magni'];
        varname=get(handles.VariableSelect,'string');
        evalin('base', [varname '=permute(' varname ',[2 1 3 4]);']);
        handles.Mi=evalin('base', varname);


        ss=size(handles.Mi);
        Image=handles.Mi(:,:,round(ss(3)/2));
        Image=permute(Image(:,end:-1:1),[2 1]);
        handles.himage1=imshow(Image,[],'parent',handles.ImageAxes);
        
        
        %varname=[UniqueID '_phase'];
        varname=get(handles.RightImageLoad,'string');
        evalin('base', [varname '=permute(' varname ',[2 1 3 4]);']);
        handles.M2=evalin('base', varname);

        Image=handles.M2(:,:,round(ss(3)/2));
        Image=permute(Image(:,end:-1:1),[2 1]);
        handles.himage2=imshow(Image,[],'parent',handles.ImageAxes2);

    case 2
        %varname=[UniqueID '_magni'];
        varname=get(handles.VariableSelect,'string');
        evalin('base', [varname '=flipdim(' varname ',2);']);
        handles.Mi=evalin('base', varname);


        ss=size(handles.Mi);
        Image=handles.Mi(:,:,round(ss(3)/2));
        Image=permute(Image(:,end:-1:1),[2 1]);
        handles.himage1=imshow(Image,[],'parent',handles.ImageAxes);

         %varname=[UniqueID '_phase'];
        varname=get(handles.RightImageLoad,'string');
        evalin('base', [varname '=flipdim(' varname ',2);']);
        handles.M2=evalin('base', varname);

        Image=handles.M2(:,:,round(ss(3)/2));
        Image=permute(Image(:,end:-1:1),[2 1]);
        handles.himage2=imshow(Image,[],'parent',handles.ImageAxes2);

    case 3

        %varname=[UniqueID '_phase'];
         varname=get(handles.RightImageLoad,'string');
       PhaseImage=double(evalin('base', varname));
        
        PhaseRange= max(PhaseImage(:))-min(PhaseImage(:));
        
        if PhaseRange>200
            PhaseRange=PhaseRange+1;
            PhaseImage=PhaseImage/PhaseRange*2*pi;

            assignin('base',varname,PhaseImage);
            handles.M2=PhaseImage;
            ss=size(handles.M2);
            Image=handles.M2(:,:,round(ss(3)/2));
            Image=permute(Image(:,end:-1:1),[2 1]);
            handles.himage2=imshow(Image,[],'parent',handles.ImageAxes2);
            set(handles.status,'string','Image scaled to 2PI.')
        else
            set(handles.status,'string','Already in the scale of 2PI.')
           return; 
        end

    case 4
        %varname=[UniqueID '_magni'];
        varname=get(handles.VariableSelect,'string');
        evalin('base', [varname '=flipdim(' varname ',1);']);
        handles.Mi=evalin('base', varname);


        ss=size(handles.Mi);
        Image=handles.Mi(:,:,round(ss(3)/2));
        Image=permute(Image(:,end:-1:1),[2 1]);
        handles.himage1=imshow(Image,[],'parent',handles.ImageAxes);

         %varname=[UniqueID '_phase'];
        varname=get(handles.RightImageLoad,'string');
        evalin('base', [varname '=flipdim(' varname ',1);']);
        handles.M2=evalin('base', varname);

        Image=handles.M2(:,:,round(ss(3)/2));
        Image=permute(Image(:,end:-1:1),[2 1]);
        handles.himage2=imshow(Image,[],'parent',handles.ImageAxes2);
end

handles.ImageAxisNum=3;
guidata(hObject, handles);



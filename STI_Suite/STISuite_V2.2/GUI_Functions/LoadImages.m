function LoadImages(hObject,handles,control)
% Wei Li, 06/14/2014
[file,file_path] = uigetfile('*.mat','MultiSelect','on');
if file==0
    return
end
my_var=load([file_path file]);
my_name=fieldnames(my_var);
warning off
switch control
    case 1
        handles.Mi=single(getfield(my_var,my_name{1}));
        ss=size(handles.Mi);
        set(handles.FrameSelect, 'Max', ss(3),'SliderStep',[1/ss(3) 2/ss(3)]); 
        Image=abs(handles.Mi(:,:,round(ss(3)/2)));
        Image=permute(Image(:,end:-1:1),[2 1]);
        handles.himage1=imshow(Image,[],'parent',handles.ImageAxes);
        set(handles.FrameSelect, 'value', round(ss(3)/2)); 
        if length(ss)==4
           set(handles.Coilslider, 'Min',1,'Max', ss(4),'SliderStep',[1/(ss(4)-1) 2/(ss(4)-1)]);
        end
        set(handles.Coilslider, 'value',1); 
        handles.ImageSize=ss;
        UniqueID=get(handles.UniqueID,'string');
        assignin('base',[UniqueID '_magni'],handles.Mi);
        handles.ImageAxisNum=3;
        set(handles.VariableSelect,'string',[UniqueID '_magni']); 
        AA=get(handles.QSMProcessing,'position');
        x0=get(handles.PhaseProcessing,'position');
        set(handles.QSMProcessing,'position',[x0(1) AA(2) AA(3) AA(4)]);
        AA=get(handles.VisualizeSelect,'position');
        set(handles.VisualizeSelect,'position',[x0(1) AA(2) AA(3) AA(4)]);
    case 2
        handles.M2=single(getfield(my_var,my_name{1}));
        ss=size(handles.M2);
        Image=handles.M2(:,:,round(ss(3)/2));
        Image=permute(Image(:,end:-1:1),[2 1]);
       
        handles.himage2=imshow(Image,[],'parent',handles.ImageAxes2);
        UniqueID=get(handles.UniqueID,'string');
        assignin('base',[UniqueID '_phase'],handles.M2);
        set(handles.RightImageLoad,'string',[UniqueID '_phase']); 
case 3
        handles.Mask=getfield(my_var,my_name{1});
        ss=size(handles.Mask);
        
        Image=abs(handles.Mask(:,:,round(ss(3)/2)));
        Image=permute(Image(:,end:-1:1),[2 1]);
        handles.himage2=imshow(Image,[],'parent',handles.ImageAxes2);
        UniqueID=get(handles.UniqueID,'string');
        assignin('base',[UniqueID '_Mask'],handles.Mask);
        set(handles.BrainMaskName,'string',[UniqueID '_Mask']); 
        set(handles.NewMaskName,'string',[UniqueID '_Mask']);
case 4
        handles.Mi=getfield(my_var,my_name{1});
        ss=size(handles.Mi);
        set(handles.FrameSelect, 'Max', ss(3),'SliderStep',[1/ss(3) 2/ss(3)]); 
        Image=abs(handles.Mi(:,:,round(ss(3)/2)));
        Image=permute(Image(:,end:-1:1),[2 1]);
        handles.himage1=imshow(Image,[],'parent',handles.ImageAxes);
        set(handles.FrameSelect, 'value', round(ss(3)/2)); 
        handles.ImageSize=ss;
case 5
        handles.M2=getfield(my_var,my_name{1});
        ss=size(handles.M2);
        Image=handles.M2(:,:,round(ss(3)/2));
        Image=permute(Image(:,end:-1:1),[2 1]);
        handles.himage2=imshow(Image,[],'parent',handles.ImageAxes2);
    case 6
        handles.Mi=getfield(my_var,my_name{1});
        ss=size(handles.Mi);
        set(handles.FrameSelect, 'Max', ss(3),'SliderStep',[1/ss(3) 2/ss(3)]); 
        Image=abs(handles.Mi(:,:,round(ss(3)/2)));
        Image=permute(Image(:,end:-1:1),[2 1]);
        handles.himage1=imshow(Image,[],'parent',handles.ImageAxes);
        set(handles.FrameSelect, 'value', round(ss(3)/2)); 
        set(handles.Coilslider, 'Min',1,'Max', ss(4),'SliderStep',[1/(ss(4)-1) 2/(ss(4)-1)]); 
        set(handles.Coilslider, 'value',1); 
        handles.ImageSize=ss;
        UniqueID=get(handles.UniqueID,'string');
        handles.ImageAxisNum=3;
        AA=get(handles.QSMProcessing,'position');
        x0=get(handles.PhaseProcessing,'position');
        set(handles.QSMProcessing,'position',[x0(1) AA(2) AA(3) AA(4)]);
        AA=get(handles.VisualizeSelect,'position');
        set(handles.VisualizeSelect,'position',[x0(1) AA(2) AA(3) AA(4)]);
        set(handles.LoadRealData,'string','***'); 
    case 7
        Imagdata=getfield(my_var,my_name{1});
        Complexdata=handles.Mi+1i*Imagdata;
        Magdata=abs(Complexdata);
        Phasedata=angle(Complexdata);
        handles.M2=Phasedata;
        handles.Mi=Magdata;
        ss=size(handles.M2);
        Image=handles.M2(:,:,round(ss(3)/2));
        Image=permute(Image(:,end:-1:1),[2 1]);
        handles.himage2=imshow(Image,[],'parent',handles.ImageAxes2);
        Image=handles.Mi(:,:,round(ss(3)/2));
        Image=permute(Image(:,end:-1:1),[2 1]);
        handles.himage1=imshow(Image,[],'parent',handles.ImageAxes);
        UniqueID=get(handles.UniqueID,'string');
        assignin('base',[UniqueID '_magni'],Magdata);
        set(handles.VariableSelect,'string',[UniqueID '_magni']); 
        assignin('base',[UniqueID '_phase'],Phasedata);
        set(handles.RightImageLoad,'string',[UniqueID '_phase']); 
        set(handles.LoadImagData,'string','***'); 
        
    case 8

        spatialres=getfield(my_var,'voxelsize');
        str1=sprintf('[ %g  %g  %g ]',spatialres);
        set(handles.LeftResolution,'string',str1); 
        set(handles.VoxelSizeForQSM,'string',str1); 
        set(handles.Voxelsize5,'string',str1); 

         TE1=getfield(my_var,'TE1');
        str1=sprintf(' %g ',TE1);
        set(handles.TE1,'string',str1); 
      
        EchoSpacing=getfield(my_var,'delta_TE');
        str1=sprintf(' %g ',EchoSpacing);
        set(handles.DeltaTE,'string',str1); 
        
        B0_vector=getfield(my_var,'B0_vector');
        str1=sprintf('[ %g  %g  %g ]',B0_vector);
        set(handles.H_Vector,'string',str1); 
        

        nEchoes=getfield(my_var,'nEchoes');
        TEs=sum(TE1+EchoSpacing*(0:nEchoes-1));
        str1=sprintf(' %g ',TEs);
        set(handles.TE_value,'string',str1); 
        
        B0Value=getfield(my_var,'B0');
        str1=sprintf(' %g ',B0Value);
        set(handles.B0Value,'string',str1); 

        
end
handles.ImageAxisNum=3;
guidata(hObject, handles);



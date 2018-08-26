function FastQSMCallBack(hObject,handles,~)

disp('Preforming Fast QSM ...')
set(handles.status,'string','Preforming Fast QSM, please wait ...')
drawnow
tic
myvarname=get(handles.ProcessedPhaseName,'string');
TissuePhase=evalin('base', myvarname);
myvarname2=get(handles.NewMaskName,'string');
NewMask=evalin('base', myvarname2);
eval(['SpatialRes=' get(handles.VoxelSizeForQSM,'string') ';'])
padsize=[8 8 20];
eval(['H=' get(handles.H_Vector,'string') ';'])
eval(['niter=' get(handles.Niter,'string') ';'])
eval(['B0=' get(handles.B0Value,'string') ';'])
eval(['TE=' get(handles.TE_value,'string') ';'])
[X0]=FastQSM(TissuePhase,NewMask,'TE',TE,'B0',B0,'H',H,'padsize',padsize,'voxelsize',SpatialRes);
UniqueID=get(handles.UniqueID,'string');
assignin('base',[UniqueID '_QSM_Fast'],X0);
handles.M2=X0;
set(handles.status,'string','Fast QSM is done.')
toc
disp('------------------------------------------------')
drawnow
guidata(hObject, handles);
SliceSelection(hObject,handles);
AdjustIntensity(hObject,handles,3)

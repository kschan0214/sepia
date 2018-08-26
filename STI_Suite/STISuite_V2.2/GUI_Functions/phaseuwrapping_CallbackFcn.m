function phaseuwrapping_CallbackFcn(hObject,handles,ctrl)

UniqueID=get(handles.UniqueID,'string');
disp('Starting unwrapping ...')
set(handles.status,'string','Starting unwrapping, please wait ...')
drawnow
phase=evalin('base', [UniqueID '_phase']);
eval(['voxelsize=' get(handles.LeftResolution,'string') ';'])
[phi Laplacian]=MRPhaseUnwrap(phase,'voxelsize',voxelsize);

UniqueID=get(handles.UniqueID,'string');
assignin('base',[UniqueID '_UnwrappedPhase'],single(phi));
assignin('base',[UniqueID '_Laplacian'],Laplacian);

set(handles.unwrappedphase,'string',[UniqueID '_UnwrappedPhase']);
handles.M2=phi;
disp('Unwrapping done.')
set(handles.status,'string','Unwrapping Done.')
guidata(hObject, handles);
SliceSelection(hObject,handles);
AdjustIntensity(hObject,handles,3)
disp('------------------------------------------------')
drawnow
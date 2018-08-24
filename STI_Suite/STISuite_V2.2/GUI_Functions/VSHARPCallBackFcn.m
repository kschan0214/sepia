function VSHARPCallBackFcn(hObject,handles,ctrl)

UniqueID=get(handles.UniqueID,'string');
set(handles.status,'string','Starting phase unwrapping and V-SHARP processing, please wait ...')
drawnow
try
    Unwrapped_Phase=evalin('base', [UniqueID '_UnwrappedPhase']);
catch
    phaseuwrapping_CallbackFcn(hObject,handles,'Unwrapping')
    Unwrapped_Phase=evalin('base', [UniqueID '_UnwrappedPhase']);
end

myvarname2=get(handles.BrainMaskName,'string');
BrainMask=evalin('base', myvarname2);

eval(['SpatialRes=' get(handles.LeftResolution,'string') ';'])

eval(['smvsize=' get(handles.PhaseIter,'string') ';'])

Unwrapped_Phase=sum(Unwrapped_Phase,4);

% [BrainMask]=SMVFiltering2(BrainMask,1.5,SpatialRes);
% BrainMask=BrainMask>0.1;

[TissuePhase,NewMask]=V_SHARP(Unwrapped_Phase,BrainMask,'voxelsize',SpatialRes,'smvsize',smvsize);

NewMask=PolishMask(NewMask);
TissuePhase=TissuePhase.*NewMask;
UniqueID=get(handles.UniqueID,'string');
assignin('base',[UniqueID '_TissuePhase_V'],TissuePhase);
assignin('base',[UniqueID '_NewMask_V'],NewMask);

set(handles.ProcessedPhaseName,'string',[UniqueID '_TissuePhase_V']);
set(handles.NewMaskName,'string',[UniqueID '_NewMask_V']);
handles.M2=TissuePhase;
set(handles.status,'string','V-SHARP Done.')
disp('------------------------------------------------')
guidata(hObject, handles);
SliceSelection(hObject,handles);
AdjustIntensity(hObject,handles,3)

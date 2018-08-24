function QSM_LSQR_Callback_Fcn(hObject,handles,ctrl)

UniqueID=get(handles.UniqueID,'string');
set(handles.status,'string','Starting QSM processing, please wait ...')
 drawnow
myvarname=get(handles.ProcessedPhaseName,'string');
TissuePhase=evalin('base', myvarname);

myvarname2=get(handles.BrainMaskName,'string');
NewMask=evalin('base', myvarname2);

eval(['SpatialRes=' get(handles.VoxelSizeForQSM,'string') ';'])
eval(['H=' get(handles.H_Vector,'string') ';'])
eval(['niter=' get(handles.Niter,'string') ';'])
eval(['B0=' get(handles.B0Value,'string') ';'])
eval(['TE=' get(handles.TE_value,'string') ';'])
eval(['Kthreshold=' get(handles.fraction_kspace,'string') ';'])
eval(['tol_step1=' get(handles.tol_step1,'string') ';'])
eval(['padsize=' get(handles.ZeroPaddingAroundBrain,'string') ';'])


params.H=H;
params.voxelsize=SpatialRes;
params.niter=niter;
params.TE=TE;
params.B0=B0;
params.tol_step1=tol_step1;
params.tol_step2=0.001;
params.Kthreshold=Kthreshold;
params.padsize=padsize;

SF=ScalingFactor(B0,TE);

[Susceptibility] = QSM_iLSQR(TissuePhase,NewMask,'params',params);
FrequencyShift=TissuePhase*SF.Freq;
assignin('base',[UniqueID '_QSM_iLSQR'],Susceptibility);
assignin('base',[UniqueID '_FrequencyShift'],FrequencyShift);

handles.Mi=FrequencyShift;
handles.M2=Susceptibility;
set(handles.status,'string','QSM done. All variables are saved in workspace, Starting with Unique ID. Have a nice day!')
disp('------------------------------------------------')
set(handles.Coilslider, 'value', 1); 
set(handles.CoilNum, 'string', '1'); 
guidata(hObject, handles);
SliceSelection(hObject,handles);
AdjustIntensity(hObject,handles,3)

function QSM_Functions(hObject,handles,ctrl)

UniqueID=get(handles.UniqueID,'string');
switch ctrl
case 'EstimateX0'
        disp('Estimate X0 ...')
        set(handles.status,'string','Quite Estimate X0, please wait ...')
         drawnow
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
        [~,X0]=EstimateX0(TissuePhase,NewMask,'TE',TE,'B0',B0,'H',H,'padsize',padsize,'voxelsize',SpatialRes);
        UniqueID=get(handles.UniqueID,'string');
        assignin('base',[UniqueID '_Chi_estimate'],X0);
        handles.M2=X0;
        set(handles.status,'string','X0 is estimated...')
        drawnow
        
    case 'QSM'
        disp('Starting QSM using iLSQR ...')
        set(handles.status,'string','Starting QSM processing, please wait ...')
         drawnow
        myvarname=get(handles.ProcessedPhaseName,'string');
        TissuePhase=evalin('base', myvarname);

        myvarname2=get(handles.NewMaskName,'string');
        NewMask=evalin('base', myvarname2);

        eval(['SpatialRes=' get(handles.VoxelSizeForQSM,'string') ';'])
 %       eval(['padsize=' get(handles.PadSizeForQSM,'string') ';'])
        eval(['H=' get(handles.H_Vector,'string') ';'])
        eval(['niter=' get(handles.Niter,'string') ';'])
        eval(['B0=' get(handles.B0Value,'string') ';'])
        eval(['TE=' get(handles.TE_value,'string') ';'])
        eval(['Kthreshold=' get(handles.fraction_kspace,'string') ';'])
    %    eval(['cropsize=' get(handles.cropsize,'string') ';'])
        eval(['tol_step1=' get(handles.tol_step1,'string') ';'])
        eval(['tol_step2=' get(handles.tol_step2,'string') ';'])

        params.H=H;
        params.voxelsize=SpatialRes;
        params.niter=niter;
        params.TE=TE;
        params.B0=B0;
        params.tol_step1=tol_step1;
        params.tol_step2=tol_step2;
        params.Kthreshold=Kthreshold;
        SF=ScalingFactor(B0,TE);
        [Susceptibility  X_LSQR] = QSM_iLSQR_Fast(TissuePhase,NewMask,'params',params);
        FrequencyShift=TissuePhase*SF.Freq;
        UniqueID=get(handles.UniqueID,'string');
        assignin('base',[UniqueID '_Chi_iLSQR'],Susceptibility);
        assignin('base',[UniqueID '_Chi_LSQR'],X_LSQR);
        assignin('base',[UniqueID '_FrequencyShift'],FrequencyShift);
        assignin('base',[UniqueID '_ScallingFactor'],SF);
        handles.Mi=FrequencyShift;
        handles.M2=Susceptibility;
        set(handles.status,'string','QSM done. All variables are saved in workspace, Starting with Unique ID. Have a nice day!')
        drawnow

case 'QSM_KIWI'
        disp(' Starting QSM using iCS ...')
        set(handles.status,'string','Starting QSM processing using KIWI, please wait ...')
        drawnow
        myvarname=get(handles.ProcessedPhaseName,'string');
        TissuePhase=evalin('base', myvarname);

        myvarname2=get(handles.NewMaskName,'string');
        NewMask=evalin('base', myvarname2);

        eval(['SpatialRes=' get(handles.VoxelSizeForQSM,'string') ';'])
        eval(['padsize=' get(handles.PadSizeForQSM,'string') ';'])
        eval(['H=' get(handles.H_Vector,'string') ';'])
        eval(['cropsize=' get(handles.cropsize,'string') ';'])
        eval(['B0=' get(handles.B0Value,'string') ';'])
        eval(['TE=' get(handles.TE_value,'string') ';'])
        
        SF=ScalingFactor(B0,TE);
        Susceptibility = QSM_iCS(TissuePhase,NewMask,'H',H,'voxelsize',SpatialRes,'padsize',padsize,'TE',TE,'B0',B0,'cropsize',cropsize);
        FrequencyShift=TissuePhase*SF.Freq;

        UniqueID=get(handles.UniqueID,'string');
        assignin('base',[UniqueID '_Susceptibility_KIWI'],Susceptibility);
        assignin('base',[UniqueID '_FrequencyShift'],FrequencyShift);
        handles.Mi=FrequencyShift;
        handles.M2=Susceptibility;
        set(handles.status,'string','QSM done. All variables are saved in workspace, Starting with Unique ID. Have a nice day!')
        drawnow
        
    case 'VSHARP'
        disp(' Starting Laplacian Phase Unwrapping and V-SHARP ...')
        set(handles.status,'string','Starting phase unwrapping and V-SHARP processing, please wait ...')
        drawnow
        Unwrapped_Phase=evalin('base', [UniqueID '_UnwrappedPhase']);

        myvarname2=get(handles.BrainMaskName,'string');
        BrainMask=evalin('base', myvarname2);

        eval(['SpatialRes=' get(handles.LeftResolution,'string') ';'])
        [TissuePhase,NewMask]=V_SHARP_Fast(Unwrapped_Phase,BrainMask,'voxelsize',SpatialRes);
        
        NewMask=PolishMask(NewMask);
        TissuePhase=TissuePhase.*NewMask;
        UniqueID=get(handles.UniqueID,'string');
        assignin('base',[UniqueID '_TissuePhase_V'],TissuePhase);
        assignin('base',[UniqueID '_NewMask_V'],NewMask);

        set(handles.ProcessedPhaseName,'string',[UniqueID '_TissuePhase_V']);
        set(handles.NewMaskName,'string',[UniqueID '_NewMask_V']);
        handles.M2=TissuePhase;
        disp(' V-SHARP done.')
        set(handles.status,'string','Phase unwrapping and V-SHARP Done. Now please fill in the QSM paramteters, and click "QSM using LSQR"')
        drawnow
end

guidata(hObject, handles);
SliceSelection(hObject,handles);
AdjustIntensity(hObject,handles,3)

function QSMEssential(hObject,handles)
UniqueID=get(handles.UniqueID,'string');
disp('------------------------------------------------')
disp(UniqueID)

disp('Starting unwrapping ...')
set(handles.status,'string','Starting unwrapping, please wait ...')
eval(['voxelsize=' get(handles.LeftResolution,'string') ';'])



%{-
drawnow
phase=evalin('base', [UniqueID '_phase']);
[phi, Laplacian]=MRPhaseUnwrap(phase,'voxelsize',voxelsize);

UniqueID=get(handles.UniqueID,'string');
assignin('base',[UniqueID '_UnwrappedPhase'],single(phi));
assignin('base',[UniqueID '_Laplacian'],Laplacian);

switch get (handles.phasePanal,'value')
    case 1
        %disp('V-SHARP Processing...')
        myvarname=get(handles.BrainMaskName,'string');
        BrainMask=evalin('base', myvarname);
        Laplacian=sum(Laplacian,4);
        NP=size(Laplacian);
        TissuePhase=zeros(NP);
        for iii=1:1
        [TissuePhase(:,:,:,iii)]=iHARPERELLA(Laplacian,BrainMask,'voxelsize',voxelsize,'niter',100);        
        set(handles.status,'string',['Echo: ' num2str(iii) ])
        drawnow;
        end
        TissuePhase=single(TissuePhase);
        assignin('base',[UniqueID '_TissuePhase'],TissuePhase);
        set(handles.status,'string','HARPERELLA Done. Now please fill in the QSM paramteters, and click "QSM using LSQR"')
        disp('------------------------------------------------')
    case 0
        %disp('V-SHARP Processing...')
        myvarname2=get(handles.BrainMaskName,'string');
        BrainMask=evalin('base', myvarname2);
        eval(['smvsize=' get(handles.PhaseIter,'string') ';'])
        Unwrapped_Phase=single(sum(phi,4));
        %[BrainMask]=SMVFiltering2(BrainMask,1.5,voxelsize);
        %BrainMask=BrainMask>0.1;
        [TissuePhase,NewMask]=V_SHARP(Unwrapped_Phase,BrainMask,'voxelsize',voxelsize,'smvsize',smvsize);
        NewMask=PolishMask(NewMask);
        TissuePhase=TissuePhase.*NewMask;
        assignin('base',[UniqueID '_TissuePhase_V'],TissuePhase);
        assignin('base',[UniqueID '_NewMask_V'],NewMask);
        BrainMask=NewMask;
        disp('------------------------------------------------')
end



eval(['H=' get(handles.H_Vector,'string') ';'])
eval(['B0=' get(handles.B0Value,'string') ';'])
eval(['TE=' get(handles.TE_value,'string') ';'])


params.H=H;
params.TE=TE;
params.B0=B0;
params.voxelsize=voxelsize;

params.niter=100;
params.tol_step1=0.01;
params.tol_step2=0.001;
params.Kthreshold=0.1;
params.padsize=[6 6 6];

SF=ScalingFactor(B0,TE);

[Susceptibility] = QSM_iLSQR(TissuePhase,BrainMask,'params',params);
FrequencyShift=TissuePhase*SF.Freq;


assignin('base',[UniqueID '_QSM_iLSQR'],Susceptibility);
assignin('base',[UniqueID '_FrequencyShift'],FrequencyShift);
Mag=evalin('base', [UniqueID '_magni']);

% eval(['TE1=' get(handles.TE1,'string') ';'])
% eval(['EchoSpacing=' get(handles.DeltaTE,'string') ';'])
%[R2Star]=R2StarMapping(Mag,TE1,EchoSpacing);
%assignin('base',[UniqueID '_QSM_R2Star'],R2Star);

handles.Mi=FrequencyShift;
handles.M2=Susceptibility;
set(handles.status,'string','QSM done. All variables are saved in workspace. Have a nice day!')
disp('------------------------------------------------')
set(handles.Coilslider, 'value', 1); 
set(handles.CoilNum, 'string', '1'); 
guidata(hObject, handles);
SliceSelection(hObject,handles);
AdjustIntensity(hObject,handles,3)

%}
filters = {'*.nii','NII-files (*.nii)'};
seed = [UniqueID '_QSM_iLSQR.nii'];
[fn,pn,~] = uiputfile(filters, sprintf('Save Workspace Variables'), seed);

fn2=[fn(1:end-14) '_FrequencyShift' '.nii'];
fn3=[fn(1:end-14) '_NewMask' '.nii'];
fn4=[fn(1:end-14) '_NewMag' '.nii'];
if fn~=0
    ezsave_nii(Susceptibility,strcat(pn,fn),[],voxelsize);
    ezsave_nii(FrequencyShift,strcat(pn,fn2),[],voxelsize);
    ezsave_nii(single(BrainMask),strcat(pn,fn3),[],voxelsize);
    ezsave_nii(Mag,strcat(pn,fn4),[],voxelsize);
    gzip(strcat(pn,fn));  delete(strcat(pn,fn))
    gzip(strcat(pn,fn2)); delete(strcat(pn,fn2))
    gzip(strcat(pn,fn3)); delete(strcat(pn,fn3))
    gzip(strcat(pn,fn4)); delete(strcat(pn,fn4))
end

        
disp('All Files Saved. Thank you for using STI Suite :)')









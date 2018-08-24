function CalliHARPERELLA(hObject,handles,ctrl)
% Wei Li, PhD

UniqueID=get(handles.UniqueID,'string');
set(handles.status,'string','Starting HARPERELLA processing, please wait ...')
drawnow

try
    Laplacian=evalin('base', [UniqueID '_Laplacian']);
catch
    phaseuwrapping_CallbackFcn(hObject,handles,'Unwrapping')
    Laplacian=evalin('base', [UniqueID '_Laplacian']);

end

myvarname=get(handles.BrainMaskName,'string');
BrainMask=evalin('base', myvarname);
Laplacian=sum(Laplacian,4);
eval(['SpatialRes=' get(handles.LeftResolution,'string') ';'])

NP=size(Laplacian);
TissuePhase=zeros(NP);
for iii=1:1
[TissuePhase(:,:,:,iii)]=iHARPERELLA(Laplacian,BrainMask,'voxelsize',SpatialRes,'niter',100);        
set(handles.status,'string',['Echo: ' num2str(iii) ])
drawnow;
end
TissuePhase=single(TissuePhase);
assignin('base',[UniqueID '_TissuePhase'],TissuePhase);
set(handles.ProcessedPhaseName,'string',[UniqueID '_TissuePhase']);
set(handles.NewMaskName,'string',[UniqueID '_Mask']);
handles.M2=TissuePhase;
set(handles.status,'string','HARPERELLA Done. Now please fill in the QSM paramteters, and click "QSM using LSQR"')
disp('------------------------------------------------')
guidata(hObject, handles);
SliceSelection(hObject,handles);
AdjustIntensity(hObject,handles,3)

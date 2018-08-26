function CalcCombinedLaplacian(hObject,handles,control)
% Wei Li, Duke University, June 2010
set(handles.status,'string','Calculate Laplacian please wait....')
drawnow
tic
str=get(handles.UniqueID,'string');
% mag=evalin('base', [str '_magni']);
% phase=evalin('base', [str '_phase']);
%BrainMask0=evalin('base', [str '_Mask']);
mag0=evalin('base', get(handles.VariableSelect,'string'));
phase0=evalin('base', get(handles.RightImageLoad,'string'));
SS00=size(phase0);
BrainMask0=evalin('base', get(handles.BrainMaskName,'string'));
eval(['voxelsize=' get(handles.LeftResolution,'string') ';'])

eval(['voxelsize=' get(handles.LeftResolution,'string') ';'])

eval(['TE1=' (get(handles.TE1,'string')) ';'])
eval(['DeltaTE=' (get(handles.DeltaTE,'string')) ';'])
TEs=TE1+DeltaTE*(0:size(phase,4)-1);

sumTE=sum(TEs);
set(handles.TE_value,'string',num2str(sumTE));



[phi Laplacian]=MEPhaseUnwrap(phase0,'voxelsize',voxelsize,'BrainMask',BrainMask0);



LP0 = MultiEchoLaplacian(phase,'voxelsize',voxelsize);

LPCombined0=sum(LP0,4);


LPCombined=zeros(SS00(1:3));

LPCombined0=LPCombined0(padsize3D(1)+1:end-padsize3D(1),padsize3D(2)+1:end-padsize3D(2),padsize3D(3)+1:end-padsize3D(3));
LPCombined(BDBX(1,1):BDBX(1,2),BDBX(2,1):BDBX(2,2),BDBX(3,1):BDBX(3,2))=LPCombined0;

LP=zeros(SS00);
LP0=LP0(padsize3D(1)+1:end-padsize3D(1),padsize3D(2)+1:end-padsize3D(2),padsize3D(3)+1:end-padsize3D(3),:);
LP(BDBX(1,1):BDBX(1,2),BDBX(2,1):BDBX(2,2),BDBX(3,1):BDBX(3,2),:)=LP0;



assignin('base',[str '_LP'],LP);
assignin('base',[str '_LPCombined'],LPCombined);
handles.M2=LPCombined;
set(handles.status,'string','Done.')
guidata(hObject, handles);
SliceSelection(hObject,handles);
AdjustIntensity(hObject,handles,3)
disp('Laplacian calculated and combined.')
toc
disp('------------------------------------------------')

function mythresh(hObject,handles)

tic

SS=size(handles.Mi);
voxelsize=[0.9 0.9 2];
try
TE1=str2double(get(handles.TE1,'string'));
DeltaTE=str2double(get(handles.DeltaTE,'string'));
catch
    TE1=1;
    DeltaTE=1;
end

try
SS(4);
catch
    SS(4)=1;
end

TEs=TE1+DeltaTE*(0:SS(4)-1);

D_30=abs(TEs-30);
nn=find(D_30==min(D_30));

Mag0=handles.Mi;
Mag=Mag0(:,:,:,nn);

SS=size(Mag);
if mod(SS(3),2)==0
    Mag(:,:,SS(3)+2)=0;
else
    Mag(:,:,SS(3)+1)=0;
end

image2=double(Mag(:,:,round(SS(3)/2)));
scaling=1/max(image2(:))*256;
image2=image2*scaling;
Mag=Mag*scaling;
try
    level=thresh_tool(permute(image2(:,end:-1:1),[2 1]),'gray');
catch
   return; 
end
set(handles.status,'string','Create masking....')
drawnow;

if isempty(level);
    return
end
NP=size(Mag);

if NP(1)>400
    Mag=Mag(1:2:end,1:2:end,:);
end
ImageMask=Mag>level;
ImageLabel = bwlabeln(ImageMask);
stats = regionprops(ImageLabel,'Area');
RegionArea = [stats.Area];
biggest = find(RegionArea==max(RegionArea));
NewMask=(ImageLabel==biggest);
STATS = regionprops(NewMask, 'FilledImage','BoundingBox');
FilledImage=STATS(1,1).FilledImage;
B=STATS(1,1).BoundingBox;
FinalMask=NewMask;
FinalMask((B(2)+0.5):(B(2)+B(5)-0.5),(B(1)+0.5):(B(1)+B(4)-0.5),(B(3)+0.5):(B(3)+B(6)-0.5))=FilledImage;
[TempMask]=SMVFiltering2(FinalMask,2,voxelsize);
Temp1=PolishMask(TempMask>0.99);
[TempMask2]=SMVFiltering2(Temp1,3,voxelsize);
TempMask2=(TempMask2>0.02);

[TempMask2]=SMVFiltering2(TempMask2,1,voxelsize);
TempMask2=TempMask2>0.97;

FinalMask0=TempMask2;

SS1=size(FinalMask);
if NP(1)>400
    FinalMask1=zeros(SS1(1)*2,SS1(2),SS1(3));
    FinalMask1(1:2:end,:,:)=FinalMask0;
    FinalMask1(2:2:end,:,:)=FinalMask0;

    FinalMask=zeros(SS1(1)*2,SS1(2)*2,SS1(3));
    FinalMask(:,1:2:end,:,:)=FinalMask1;
    FinalMask(:,2:2:end,:,:)=FinalMask1;
    
    FinalMask=FinalMask(1:NP(1),1:NP(2),:);
end


if mod(SS(3),2)==0
    FinalMask=FinalMask(:,:,1:end-2);
else
    FinalMask=FinalMask(:,:,1:end-1);
end

 BrainMask=FinalMask;
 
 
set(handles.status,'string',['Echo: ' num2str(nn) ])
drawnow;

set(handles.status,'string','Done.')
FinalMask=single(BrainMask);
handles.M2=FinalMask;
handles.FinalMask=FinalMask;
try
UniqueID=get(handles.UniqueID,'string');
catch
UniqueID='S001';
end
    
assignin('base',[UniqueID '_Mask'],single(FinalMask));
try
set(handles.BrainMaskName,'string',[UniqueID '_Mask']); 
set(handles.NewMaskName,'string',[UniqueID '_Mask']);
end

ss=size(handles.M2);

Image=abs(handles.M2(:,:,round(ss(3)/2)));
Image=permute(Image(:,end:-1:1),[2 1]);

handles.himage2=imshow(Image,[],'parent',handles.ImageAxes2);
set(handles.Coilslider, 'value', 1); 
set(handles.CoilNum, 'string', '1'); 
toc
disp('------------------------------------------------')
guidata(hObject, handles);




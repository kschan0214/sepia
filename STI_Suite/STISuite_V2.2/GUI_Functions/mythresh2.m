function mythresh(hObject,handles)

try
    level=thresh_tool(get(handles.himage1,'CData'),'gray');
catch
   return; 
end

if isempty(level);
    return
end
Mag=handles.Mi(:,:,:,get(handles.Coilslider,'value'));

SS=size(Mag);

if SS(1)>400 
    Mag=Mag(1:2:end,1:2:end,:,4);
end

tic
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
[TempMask]=SMVFiltering(FinalMask,2);
Temp1=PolishMask(TempMask>0.99);
[TempMask2]=SMVFiltering(Temp1,3);
TempMask2=(TempMask2>0.02);

[TempMask2]=SMVFiltering(TempMask2,1);
TempMask2=TempMask2>0.97;

FinalMask=zeros(SS(1),size(FinalMask,2),size(FinalMask,3))=TempMask2;

FinalMask=zeros(SS(1),size(FinalMask,2),size(FinalMask,3))=TempMask2;
FinalMask=TempMask2;

toc
handles.M2=FinalMask;
handles.FinalMask=FinalMask;
UniqueID=get(handles.UniqueID,'string');
assignin('base',[UniqueID '_Mask'],double(FinalMask));

set(handles.BrainMaskName,'string',[UniqueID '_Mask']); 
set(handles.NewMaskName,'string',[UniqueID '_Mask']);
ss=size(handles.M2);

Image=abs(handles.M2(:,:,round(ss(3)/2)));
Image=permute(Image(:,end:-1:1),[2 1]);

handles.himage2=imshow(Image,[],'parent',handles.ImageAxes2);
set(handles.Coilslider, 'value', 1); 
set(handles.CoilNum, 'string', '1'); 
disp('------------------------------------------------')
guidata(hObject, handles);




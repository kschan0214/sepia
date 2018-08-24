
% voxelsize=[0.9 0.9 2];
% ImageMask=Mag>level;
ImageLabel = bwlabeln(temp1all(:,:,:,6));
stats = regionprops(ImageLabel,'Area');
RegionArea = [stats.Area];
biggest = find(RegionArea==max(RegionArea));
NewMask=(ImageLabel==biggest);
STATS = regionprops(NewMask, 'FilledImage','BoundingBox');
FilledImage=STATS(1,1).FilledImage;
B=STATS(1,1).BoundingBox;
FinalMask=NewMask;

% STATS = regionprops(NewMask, 'FilledImage','BoundingBox');
% FilledImage=STATS(1,1).FilledImage;
% B=STATS(1,1).BoundingBox;
% FinalMask=NewMask;
% FinalMask((B(2)+0.5):(B(2)+B(5)-0.5),(B(1)+0.5):(B(1)+B(4)-0.5),(B(3)+0.5):(B(3)+B(6)-0.5))=FilledImage;
% [TempMask]=SMVFiltering2(FinalMask,2,voxelsize);
% Temp1=PolishMask(TempMask>0.99);
% [TempMask2]=SMVFiltering2(Temp1,3,voxelsize);
% TempMask2=(TempMask2>0.02);
% 
% [TempMask2]=SMVFiltering2(TempMask2,1,voxelsize);
% TempMask2=TempMask2>0.97;
% 
% FinalMask0=TempMask2;
% 
% SS1=size(FinalMask);
% if NP(1)>400
%     FinalMask1=zeros(SS1(1)*2,SS1(2),SS1(3));
%     FinalMask1(1:2:end,:,:)=FinalMask0;
%     FinalMask1(2:2:end,:,:)=FinalMask0;
% 
%     FinalMask=zeros(SS1(1)*2,SS1(2)*2,SS1(3));
%     FinalMask(:,1:2:end,:,:)=FinalMask1;
%     FinalMask(:,2:2:end,:,:)=FinalMask1;
%     
%     FinalMask=FinalMask(1:NP(1),1:NP(2),:);
% end
% 
% 
% if mod(SS(3),2)==0
%     FinalMask=FinalMask(:,:,1:end-2);
% else
%     FinalMask=FinalMask(:,:,1:end-1);
% end
% 
%  BrainMask(:,:,:,nn)=FinalMask;
% set(handles.status,'string',['Echo: ' num2str(nn) ])
% drawnow;
% for iii=nn-1:-1:1
%     phi0_Low=SMVFiltering2(BrainMask(:,:,:,iii+1),2,voxelsize);
%     phi0_Low=phi0_Low>0.1;
%     Mag=Mag0(:,:,:,iii).*phi0_Low;
%     
%     SS=size(Mag);
%     if mod(SS(3),2)==0
%         Mag(:,:,SS(3)+2)=0;
%     else
%         Mag(:,:,SS(3)+1)=0;
%     end
%     NP=size(Mag);
%     if NP(1)>400
%         Mag=Mag(1:2:end,1:2:end,:);
%     end
%     ImageMask=Mag>levels(iii);
%     ImageLabel = bwlabeln(ImageMask);
%     stats = regionprops(ImageLabel,'Area');
%     RegionArea = [stats.Area];
%     biggest = find(RegionArea==max(RegionArea));
%     NewMask=(ImageLabel==biggest);
%     STATS = regionprops(NewMask, 'FilledImage','BoundingBox');
%     FilledImage=STATS(1,1).FilledImage;
%     B=STATS(1,1).BoundingBox;
%     FinalMask=NewMask;
%     FinalMask((B(2)+0.5):(B(2)+B(5)-0.5),(B(1)+0.5):(B(1)+B(4)-0.5),(B(3)+0.5):(B(3)+B(6)-0.5))=FilledImage;
%     [TempMask]=SMVFiltering2(FinalMask,2,voxelsize);
%     Temp1=PolishMask(TempMask>0.99);
%     [TempMask2]=SMVFiltering2(Temp1,3,voxelsize);
%     TempMask2=(TempMask2>0.02);
% 
%     [TempMask2]=SMVFiltering2(TempMask2,1,voxelsize);
%     TempMask2=TempMask2>0.97;
% 
%     FinalMask0=TempMask2;
% 
%     tic
%     if NP(1)>400
%         FinalMask1=zeros(SS1(1)*2,SS1(2),SS1(3));
%         FinalMask1(1:2:end,:,:)=FinalMask0;
%         FinalMask1(2:2:end,:,:)=FinalMask0;
% 
%         FinalMask=zeros(SS1(1)*2,SS1(2)*2,SS1(3));
%         FinalMask(:,1:2:end,:,:)=FinalMask1;
%         FinalMask(:,2:2:end,:,:)=FinalMask1;
% 
%         FinalMask=FinalMask(1:NP(1),1:NP(2),:);
%     end
% 
%     if mod(SS(3),2)==0
%         FinalMask=FinalMask(:,:,1:end-2);
%     else
%         FinalMask=FinalMask(:,:,1:end-1);
%     end
% 
%     BrainMask(:,:,:,iii)=FinalMask;
%     
% set(handles.status,'string',['Echo: ' num2str(iii) ])
% drawnow;
% end
% set(handles.status,'string','Done.')
% FinalMask=BrainMask;
% handles.M2=FinalMask;
% handles.FinalMask=FinalMask;
% UniqueID=get(handles.UniqueID,'string');
% assignin('base',[UniqueID '_Mask'],double(FinalMask));
% 
% set(handles.BrainMaskName,'string',[UniqueID '_Mask']); 
% set(handles.NewMaskName,'string',[UniqueID '_Mask']);
% ss=size(handles.M2);
% 
% Image=abs(handles.M2(:,:,round(ss(3)/2)));
% Image=permute(Image(:,end:-1:1),[2 1]);
% 
% handles.himage2=imshow(Image,[],'parent',handles.ImageAxes2);
% set(handles.Coilslider, 'value', 1); 
% set(handles.CoilNum, 'string', '1'); 
% toc
% disp('------------------------------------------------')
% guidata(hObject, handles);




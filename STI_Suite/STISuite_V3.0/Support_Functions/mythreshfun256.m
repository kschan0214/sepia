function  [BrainMask1 DataIndex]=mythreshfun256(img,image,DataIndex0)

if nargin==1
img=padarray(img,[0 18 0]);
img(:,1:18,:)=img(:,end-17:end,:);
img(:,end-17:end,:)=img(:,1:18,:);

img1=sum(sum(img,1),2);
z_index= (img1==max(img1(:)));
imagez=img(:,:,z_index);
imagez=imagez(:);
meanthres=mean(img(:))*0.6;
[n,xout]=hist(imagez,100);
nout=(xout>meanthres).*n;
MaskThreshold=xout(nout==max(nout))*0.6;
ImageMask=img>MaskThreshold;

ImageLabel = bwlabeln(ImageMask);
stats = regionprops(ImageLabel,'Area');
RegionArea = [stats.Area];
biggest = find(RegionArea==max(RegionArea));
NewMask=(ImageLabel==biggest);
STATS = regionprops(NewMask, 'FilledImage','BoundingBox');
FilledImage=STATS(1,1).FilledImage;
B=STATS(1,1).BoundingBox;
BrainMask=NewMask;
BrainMask((B(2)+0.5):(B(2)+B(5)-0.5),(B(1)+0.5):(B(1)+B(4)-0.5),(B(3)+0.5):(B(3)+B(6)-0.5))=FilledImage;

STR1=ball3D(3);
STR1(4,4,4)=1;
BrainMask1=imerode(BrainMask,STR1);

ImageLabel = bwlabeln(BrainMask1);
stats = regionprops(ImageLabel,'Area');
RegionArea = [stats.Area];
biggest = find(RegionArea==max(RegionArea));
NewMask=(ImageLabel==biggest);
STATS = regionprops(NewMask, 'FilledImage','BoundingBox');
FilledImage=STATS(1,1).FilledImage;
B=STATS(1,1).BoundingBox;
FinalMask=NewMask;
FinalMask((B(2)+0.5):(B(2)+B(5)-0.5),(B(1)+0.5):(B(1)+B(4)-0.5),(B(3)+0.5):(B(3)+B(6)-0.5))=FilledImage;
BrainMask1=imdilate(FinalMask,STR1);

Mask00=sum(sum(BrainMask1,2),3);
Index=find(Mask00>0);
Datamiddle=round(mean(Index));
DataIndex=[Datamiddle-108 Datamiddle+109];
if DataIndex(1)<1
    DataIndex=[1 218];
elseif DataIndex(end)>256
     DataIndex=[39 256];
end

BrainMask1=BrainMask1(DataIndex(1):DataIndex(2),:,:);

end

if nargin==3
    image=padarray(image,[0 18 0]);
    BrainMask1=image(DataIndex0(1):DataIndex0(2),:,:);
end

    




% B3=ball3D(1);
% B3(2,2,2)=1;
% BrainMask1=imdilate(BrainMask1,B3);

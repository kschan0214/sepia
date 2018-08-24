function  BrainMask1=mythreshfun(img)

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


% B3=ball3D(1);
% B3(2,2,2)=1;
% BrainMask1=imdilate(BrainMask1,B3);

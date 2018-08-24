function FinalMask = CreateMask(Mag,level)

Mag = imscale(abs(Mag));
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

end

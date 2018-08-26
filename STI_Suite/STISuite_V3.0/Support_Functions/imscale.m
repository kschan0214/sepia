function imsc = imscale(im)

maxv = max((im(:)));
minv = min((im(:)));
imsc = ((im)-minv)./(maxv-minv);
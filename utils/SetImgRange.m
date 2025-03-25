%% uility function
function res = SetImgRange(img,range)
    imgMax = range(2);
    imgMin = range(1);
    img(img<imgMin) = imgMin;
    img(img>imgMax) = imgMax;
    img(isinf(img)) = imgMin;
    img(isnan(img)) = imgMin;
    res = img;
end
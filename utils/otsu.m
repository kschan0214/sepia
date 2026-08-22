% From:
% https://en.wikipedia.org/wiki/Otsu%27s_method
function level = otsu(histogramCounts)
total = sum(histogramCounts); % total number of pixels in the image 
%% OTSU automatic thresholding
top = length(histogramCounts);
sumB = 0;
wB = 0;
maximum = 0.0;
sum1 = dot(0:top-1, histogramCounts);
for ii = 1:top
    wF = total - wB;
    if wB > 0 && wF > 0
        mF = (sum1 - sumB) / wF;
        val = wB * wF * ((sumB / wB) - mF) * ((sumB / wB) - mF);
        if ( val >= maximum )
            level = ii;
            maximum = val;
        end
    end
    wB = wB + histogramCounts(ii);
    sumB = sumB + (ii-1) * histogramCounts(ii);
end
end
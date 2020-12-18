function [mapImg, cRange] = imageScale(cMap, cRange)
% quick command to plot image and set NaNs to be transparent. cRange
% defines a symmetric color range, based on top 'cRange' percentile in the image.
% usage: [mapImg, cRange] = imageScale(cMap, cRange)

if ~exist('cRange', 'var')
    cRange = abs(prctile(cMap(:),97.5));
end

cRange = [-cRange cRange];
mapImg = imshow(cMap,cRange);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.

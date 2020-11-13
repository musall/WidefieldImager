function [areaMask, areaLabels] = getAreaMasks(mask, minSize)
%code to return individual areas from the allen atlas. Returns logical
%masks for each area on both hemisspheres and the correpsonding label.

load('allenDorsalMapSM.mat')
allenMask = dorsalMaps.allenMask;

if ~exist('mask','var') || isempty(mask)
    mask = false(size(allenMask)); %reject pixels within this mask
end

if ~exist('minSize','var') || isempty(minSize)
    minSize = 0; %minimum area size
end

areaMask = {}; areaLabels = {};
leftIdx = find(ismember(dorsalMaps.sidesSplit,'L'));
rightIdx = find(ismember(dorsalMaps.sidesSplit,'R'));
Cnt = 0;
for iAreas = 1 : length(leftIdx)
    cIdx = [leftIdx(iAreas) rightIdx(iAreas)];
    cOutline = dorsalMaps.edgeOutlineSplit{leftIdx(iAreas)};
    cMask = poly2mask(cOutline(:,2), cOutline(:,1),size(allenMask,1),size(allenMask,2));
    cOutline = dorsalMaps.edgeOutlineSplit{rightIdx(iAreas)};
    cMask = cMask | poly2mask(cOutline(:,2), cOutline(:,1),size(allenMask,1),size(allenMask,2));
    cMask = arrayCrop(cMask,mask);
    
    if nansum(cMask(:)) > minSize
        Cnt = Cnt + 1;
        areaMask{Cnt} = cMask;
        areaLabels{Cnt} = dorsalMaps.labelsSplit{cIdx(1)};
    end
end
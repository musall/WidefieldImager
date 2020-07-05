function [Vout, regC, T, hemoVar] = SvdHemoCorrect(U, blueV, hemoV, sRate, frameCnt, highCut, smoothBlue)
% function [Vout, regC, T] = SvdHemoCorrect(U, blueV, hemoV, sRate, smoothBlue)
% 
% Does local hemodynamic correction for widefield imaging, in SVD space.
% U and V is an SVD representation of is the neural signal you are trying to correct
%
% hemoV is the other non-neural signals that you are using to measure hemodynmaics 
% It should be compressed by the same U.
%
% sRate is sampling frequency. 
%
% Outputs: 
% Vout is corrected signal
% regC are regression coefficients for each pixel to match violet to blue channel
% T is transformation matrix that predicts V from hemoV
% hemoVar is the variance in blue that is explained by the violet signal

if ~exist('highCut', 'var') || isempty(highCut)
    highCut = 10; %upper frequency threshold for low-pass smoothing filter
end

if ~exist('smoothBlue', 'var') || isempty(smoothBlue)
    smoothBlue = false; %if blue channel should be smoothed
end

%% pre-process V and U and apply mask to U
[A,B] = size(blueV);
blueV = reshape(blueV,A,[])';
hemoV = reshape(hemoV,A,[])';

% subtract means
blueV = bsxfun(@minus, blueV, nanmean(blueV));
hemoV = bsxfun(@minus, hemoV, nanmean(hemoV));

% high-pass blueV and hemoV above 0.1Hz
[b, a] = butter(2,0.1/sRate, 'high');
blueV(~isnan(blueV(:,1)),:) = single(filtfilt(b,a,double(blueV(~isnan(blueV(:,1)),:))));
hemoV(~isnan(blueV(:,1)),:) = single(filtfilt(b,a,double(hemoV(~isnan(blueV(:,1)),:))));

% get core pixels from U
mask = isnan(U(:,:,1));
U = arrayShrink(U,mask,'merge'); %only use selected pixels from mask

%% smooth hemo V
if sRate > highCut
    hemoV = smoothWidefield(hemoV,frameCnt,sRate,highCut); %smooth violet channel
    if smoothBlue
        hemoV = smoothWidefield(hemoV,frameCnt,sRate,highCut); %smooth blue channel
    end
end

%% compute single pixel time traces and regression coefficients. Always 500 at a time to prevent memory issues.
regC = zeros(1,size(U,1),'single');
ind = 0:500:size(U,1);
for x = 1:length(ind)  
    if x == length(ind)
        a = (U(ind(x)+1:end,:) * blueV');
        b = (U(ind(x)+1:end,:) * hemoV');
        regC(ind(x)+1:end) = nansum(a.*b, 2) ./ nansum(b.*b, 2);
    else
        a = (U(ind(x)+1:ind(x+1),:) * blueV');
        b = (U(ind(x)+1:ind(x+1),:) * hemoV');
        regC(ind(x)+1:ind(x+1)) = nansum(a.*b,2) ./ nansum(b.*b,2);
    end
end
clear a b

%% compute the corresponding V-space transformation matrix
T = pinv(U) * bsxfun(@times, regC(:), U);

% make the prediction
Vout = blueV - hemoV*T';
Vout = bsxfun(@minus, Vout, nanmean(Vout,1)); %subtract mean

%% compute variance explained
f1Pow = nansum(blueV(:).^2);
f1Powcor = nansum(Vout(:).^2);
hemoVar = 100*(f1Pow-f1Powcor)/f1Pow;
fprintf('%f percent variance explained by hemo signal\n', hemoVar);

% Transpose to return conventional nSVs x nTimes output
Vout = reshape(Vout', A, B);


%% nested functions
function data = smoothWidefield(data,frameCnt,sRate,highCut)
[b, a] = butter(2,highCut/sRate, 'low'); %filter below 15Hz to smooth data
for iTrials = 1:size(frameCnt,2)
    if iTrials == 1
        cIdx = 1 : frameCnt(iTrials); % first trial
    else
        cIdx = sum(frameCnt(1:iTrials-1)) + 1 : sum(frameCnt(1:iTrials)); % other trials
    end
    
    cData = data(cIdx, :); %get data for current trial
    nanIdx = ~isnan(cData(:,1)); %make sure to only use non-NaN frames
    cData = cData(nanIdx,:);
    cData = [repmat(cData(1,:),10,1); cData; repmat(cData(end,:),10,1)];
    cData = single(filtfilt(b,a,double(cData)))';
    data(cIdx(nanIdx),:) = cData(:, 11:end-10)';
end

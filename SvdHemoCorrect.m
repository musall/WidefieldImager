function [Vout, regC, T, hemoVar] = SvdHemoCorrect(U, blueV, hemoV, sRate, smoothBlue)
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

if ~exist('smoothBlue', 'var') || isempty(smoothBlue)
    smoothBlue = false;
end

%% pre-process V and U and apply mask to U
[A,B,C] = size(blueV);
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
if sRate > 15
    [b, a] = butter(2,15/sRate, 'low');
    blueV = reshape(blueV',A,B,C);
    hemoV = reshape(hemoV',A,B,C);
    
    for iTrials = 1:C
        cIdx = ~isnan(hemoV(1,:,iTrials)); %make sure to only use non-NaN frames
        if smoothBlue
            temp = blueV(:,cIdx,iTrials)';
            temp = [repmat(temp(1,:),10,1); temp; repmat(temp(end,:),10,1)];
            temp = single(filtfilt(b,a,double(temp)))';
            blueV(:,cIdx,iTrials) = temp(:, 11:end-10);
        end
        
        temp = hemoV(:,cIdx,iTrials)';
        temp = [repmat(temp(1,:),10,1); temp; repmat(temp(end,:),10,1)];
        temp = single(filtfilt(b,a,double(temp)))';
        hemoV(:,cIdx,iTrials) = temp(:, 11:end-10);
    end
    blueV = reshape(blueV,A,[])';
    hemoV = reshape(hemoV,A,[])';
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
Vout = reshape(Vout', A, B, C);

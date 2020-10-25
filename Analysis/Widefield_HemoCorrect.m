function data = Widefield_HemoCorrect(data,hemoData,baseline,smoothFact)
% This code is to perform hemodynamic correction on the raw wiedefield
% data. Mostly used to correct single trial data. 
% 
% data is the channel that contains neural+intrinsic signal, e.g. with blue
% illumination with GCaMP imaging. hemoData is the channel that contains
% only intrinsic signal. both should be 3D matrices of size X x Y x frames.
% X and Y are dimensions in pixel sace.
% baseline is a index vector to determine which frames to use for baseline
% correction. Data will be divided and subtracted by average data over all
% 'baseline' frames.
% smoothFact determines the amount of smoothing that is applied to hemoData
% before correction. Smoothing is a boxcar filter of size smoothFact^2. At
% 30Hz sampling rate, a size of 5 has yielded good results.

data = single(data); 
[A,B,C] = size(data);
dataAvg = mean(data(:,:,baseline),3);
data = bsxfun(@minus, data, dataAvg); % subtract baseline mean
data = bsxfun(@rdivide, data, dataAvg); % divide by baseline mean
data = reshape(data,[],C);

hemoData = single(hemoData);
hemoAvg = mean(hemoData(:,:,baseline),3);
hemoData = bsxfun(@minus, hemoData, hemoAvg); % subtract baseline mean
hemoData = bsxfun(@rdivide, hemoData, hemoAvg); % divide by baseline mean
hemoData = reshape(hemoData,[],C);

%% smooth hemo data
smoothFact = smoothFact-1+mod(smoothFact,2); % ensure kernel length is odd
n = size(hemoData,2);
cbegin = cumsum(hemoData(:,1:smoothFact-2),2);
cbegin = bsxfun(@rdivide, cbegin(:,1:2:end), 1:2:(smoothFact-2));
cend = cumsum(hemoData(:,n:-1:n-smoothFact+3),2);
cend = bsxfun(@rdivide, cend(:,end:-2:1), (smoothFact-2:-2:1));

hemoData = conv2(reshape(hemoData,[],C),ones(1,smoothFact)/smoothFact,'full'); %smooth trace with moving average of 'smoothFact' points
hemoData = [cbegin hemoData(:,smoothFact:end-smoothFact+1) cend];

%% perform regression and scale hemo to blue channel
if isa(data,'gpuArray')
    m = zeros(size(data,1),1,'gpuArray');
    b = zeros(size(data,1),1,'gpuArray');
else
    m = zeros(size(data,1),1);
    b = zeros(size(data,1),1);
end

warning('off','MATLAB:rankDeficientMatrix');
for iPix = 1:size(data,1)
  theta = [hemoData(iPix,:)' ones(size(data,2),1)] \ data(iPix,:)';
  m(iPix) = theta(1);
  b(iPix) = theta(2);
end
warning('on','MATLAB:rankDeficientMatrix');

hemoData = bsxfun(@times, hemoData, m); % multiply hemoData with regression coefficient
hemoData = bsxfun(@plus, hemoData, b); % add constant offset

% [corrMatrix,m] = regression(bsxfun(@minus, hemoData, mean(hemoData,2)),bsxfun(@minus, data, mean(data,2))); %run regression on 
% corrMatrix = reshape(corrMatrix,A,B); %output correlation matrix

%% subtract hemo from blue, reshape and ensure there is no baseline offset
data = data - hemoData; %subtract scaled hemoChannel from data
data = reshape(data,A,B,C);
dataAvg = mean(data(:,:,baseline),3);
data = bsxfun(@minus, data, dataAvg); % correct baseline offset


% This example code demonstrates how to analyze raw imaging data by spatial
% downsampling and baseline-correction. It shows the results by plotting a
% map of stimulus-triggered activity and a trace of fluorescence change
% over all pixels.
%
% You can test this on either behavioral or mapping example data. Just
% change the variable 'fPath' to the datapath that contains the imaging
% data. For questions, contact simon.musall@gmail.com

%% construct path to data folder and give some basic info
opts.fPath = fPath; %path to imaging data
opts.fName = 'Frames_2_640_540_uint16'; %name of imaging data files.
opts.stimLine = 4; %analog line that contains stimulus trigger.
opts.trigLine = [2 3]; %analog lines for blue and violet light triggers.
opts.preStim = 1; %pre-stimulus duration in seconds
opts.postStim = 2; %post-stimulus duration in seconds
opts.plotChans = false; %flag to show separate channels when loading data
opts.sRate = 30; %sampling rate in Hz
opts.downSample = 4; %spatial downsampling factor
opts.hemoCorrect = true; %hemodynamic correction is optional

%% load imaging data
rawFiles = dir([opts.fPath filesep opts.fName '*']); %find data files
load([opts.fPath filesep 'frameTimes_0001.mat'], 'imgSize') %get size of imaging data
imgSize = floor(imgSize ./ opts.downSample); %adjust for downsampling
stimOn = (opts.preStim*opts.sRate); %frames before stimulus onset
baselineDur = 1 : min([opts.sRate stimOn]); %use the first second or time before stimulus as baseline

nrFrames = (opts.preStim + opts.postStim) * opts.sRate; %frames per trial

nrTrials = length(rawFiles); %nr of trials
% nrTrials = 5; %make this quicker to run by using only a part of all trials

allData = NaN(imgSize(1),imgSize(2),nrFrames, nrTrials, 'single'); %pre-allocate data array for all trials
for trialNr = 1 : nrTrials
    
    [~,~,a] = fileparts(rawFiles(trialNr).name); %get data type (should be .dat for raw or .mj2 for compressed data)
    [bData,~,vData] = splitChannels(opts,(trialNr),a);
    
    if opts.hemoCorrect
        data = Widefield_HemoCorrect(bData,vData,baselineDur,5); %hemodynamic correction
    else
        data = bData; %use blue data
    end
            
    %spatially downsample imaging data
    data = arrayResize(data, 4); %this reduces resolution to ~80um / pixel
    allData(:,:,:,trialNr) = data;
    
    if rem(trialNr,floor(nrTrials / 5)) == 0
        fprintf('%d / %d files loaded\n', trialNr, nrTrials)
    end
end
clear data bData vData

if ~opts.hemoCorrect %this step is automatically done during hemodynamic correction
    
    % compute dF/F by subtracting and dividing the pre-stimulus baseline
    baselineAvg = nanmean(nanmean(allData(:,:, baselineDur, :),3), 4);
    allData = reshape(allData, imgSize(1),imgSize(2), []); %merge all frames to subtract and divide baseline
    allData = bsxfun(@minus, allData, dataAvg); % subtract baseline
    allData = bsxfun(@rdivide, allData, dataAvg); % divide baseline
    allData = reshape(allData, imgSize(1),imgSize(2),nrFrames,nrTrials); %shape back to initial form
    
end

%% show an example figure for stimulus triggered activity
figure

%show an activity map
colorRange = 0.03; %range of colorscale for dF/F

subplot(1,2,1);
avgMap = nanmean(nanmean(allData(:,:, stimOn + 1 : end, :),3),4); %show average activity after stimulus onset
imageScale(avgMap, colorRange); %show average activity after stimulus onset
colormap(colormap_blueblackred(256)); colorbar
title('Stimulus-triggered activity')

%show an activity trace
subplot(1,2,2); hold on;
meanTrace = squeeze(nanmean(reshape(allData,[], size(allData,3), size(allData,4)),1)); %average activity over all pixels
timeTrace = ((1:nrFrames) ./ opts.sRate) - opts.preStim; %time in seconds
plotLine = stdshade(meanTrace', 0.5, 'r', timeTrace); %plot average activity
plot([0 0], plotLine.Parent.YLim, '--k'); %show stimulus onset
plot(plotLine.Parent.XLim, [0 0], '--k'); %show zero line
xlim([min(timeTrace) max(timeTrace)])
axis square; xlabel('time after stimulus (s)'); ylabel('fluorescence change (dF/F)');
title('Average change over all pixels')


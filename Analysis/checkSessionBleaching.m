% checkSessionBleaching

%% Set basic variables
fclose('all');
cPath = 'G:\Mapping\Animals\Dummy Subject\Default\23-Aug-2020_2'; %Widefield data path
fName = 'Frames'; %name imaging files
numChans = 1; %number of excitation wavelengths (1 for blue only, 2 for alternating illumination)
dataType = 'uint8'; %type of imaging data usually uint16 or uint8

%% load data
rawVids = dir([cPath filesep fName '_*']); %video files

% get trial numbers for each file and sort in ascending order
trials = zeros(1,length(rawVids));
for x = 1 : length(rawVids)
    [~,a] = fileparts(rawVids(x).name);
    a = textscan(a,'%s','delimiter','_');
    trials(x) = str2double(a{1}{end});
end
[trials,sortIdx] = sort(trials,'ascend');

% get single frame data
allData = cell(1,length(trials));
allTimes = cell(1,length(trials));
for iTrials = 1 : length(trials)
    
    if numChans == 1
        load([cPath filesep 'frameTimes_' num2str(trials(iTrials))], 'imgSize', 'frameTimes') %get size of imaging file
        cFile = [cPath filesep rawVids(sortIdx(iTrials)).name];
        [~, cData] = loadRawData(cFile, 'Frames', dataType, imgSize);
        if length(size(cData)) == 4 % convert rgb image to gray
            cData =  squeeze(sum(cat(3, 0.2989 .* cData(:,:,1,:), 0.5870 .* cData(:,:,2,:), 0.1140 .* cData(:,:,3,:)), 3));
        end
    end
    allData{iTrials} = mean(reshape(cData,[], size(cData,3),1))'; %keep frame average
    allTimes{iTrials} = frameTimes;
end

% combine all trials
allData = cat(1, allData{:});
allTimes = cat(1, allTimes{:});
allTimes = (allTimes - allTimes(1)) .* 1440; %convert to minutes
allData(zscore(diff(allTimes)) > 4) = NaN; %to avoid lines between trials

%% show plot
figure('name', 'check photobleaching')
plot(allTimes,allData, 'linewidth', 2, 'Color', 'k'); 
axis square; xlabel('Time (minutes)'); ylabel('Mean fluorescence');
title(sprintf('Mean frame intensity\ndata path: %s', cPath))

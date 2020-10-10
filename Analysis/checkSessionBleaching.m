function [allData, allTimes] = checkSessionBleaching(cPath)
% this code can be used to quickly assess changes in image intensitiy due
% to photobleaching. it calculates the mean intensity in each frame and
% plots it as a function of time for the entire session.

%% Set basic variables
fclose('all');
fName = 'Frames'; %name format for imaging files
numChans = 1; %number of excitation wavelengths (1 for blue only, 2 for alternating illumination)
dataType = 'uint16'; %type of imaging data usually uint16 or uint8

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
        cFile = [cPath filesep rawVids(sortIdx(iTrials)).name];
        try
            load([cPath filesep 'frameTimes_' num2str(trials(iTrials))], 'imgSize', 'frameTimes') %get size of imaging file
            [~, cData] = loadRawData(cFile, 'Frames', dataType, imgSize);
        catch
            [header, cData] = loadRawData(cFile, 'Frames', dataType); %if file was written was old version, the image size is stored in the binary file
            if isempty(header) %if no header is found frametimes should be
                load([cPath filesep 'frameTimes_' num2str(trials(iTrials))], 'frameTimes');
            else
                frameTimes = header(1:end-4);
            end
        end
        if length(size(cData)) == 4
            if size(cData,3) == 3 % convert rgb image to gray
                cData =  squeeze(sum(cat(3, 0.2989 .* cData(:,:,1,:), 0.5870 .* cData(:,:,2,:), 0.1140 .* cData(:,:,3,:)), 3));
            else
                cData = squeeze(cData(:,:,1,:));
            end
        end
        % dont use first and last frames
        cData = cData(:,:,16:end-15);
        frameTimes = frameTimes(16:end-15);
    end
    allData{iTrials} = mean(reshape(cData,[], size(cData,3),1))'; %keep frame average
    allTimes{iTrials} = frameTimes;
end

% combine all trials
allData = cat(1, allData{:});
allTimes = cat(1, allTimes{:});
allTimes = (allTimes - allTimes(1)) .* 1440; %convert to minutes
allData(zscore(diff(allTimes)) > 1) = NaN; %to avoid lines between trials
allData(~isnan(allData)) = smooth(allData(~isnan(allData)), 50); % do some smooth

%% show plot
figure('name', 'check photobleaching')
plot(allTimes,allData, 'linewidth', 2, 'Color', 'k'); 
axis square; xlabel('Time (minutes)'); ylabel('Mean fluorescence');
title(sprintf('Mean frame intensity\ndata path: %s', cPath))

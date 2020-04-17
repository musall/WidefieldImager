function [blueData,blueTimes,hemoData,hemoTimes,stimOn,falseAlign,blueFrames,offset] = splitChannelsPart(opts,trialNr,blueFrames,offset)
% Code to separate blue and violet channel from widefield data. This needs
% analog data that contains a stimulus onset from which a pre and
% poststimulus dataset can be returned. Also requires trigger channels for
% blue/violet LEDs and some dark frames at the end.
% This code is the newer version of Widefield_CheckChannels.

% get widefield data
cFile = [opts.fPath filesep opts.fName '_' num2str(trialNr) '.dat']; %current file to be read
[header,data] = Widefield_LoadData(cFile,'Frames'); %load video data
sRate = opts.frameRate;
falseAlign = false;

if ~exist('blueFrames','var') || isempty(blueFrames) 
    %reshape data to compute mean frame intensities
    dSize = size(data);
    data = reshape(data,[],dSize(end));
    temp = zscore(median(single(data)));
    data = squeeze(reshape(data,dSize));
    frameTimes = header(1:end-length(dSize)) * (86400*1e3); %extract frame times from header and convert to millisecond timestamps
    
    %% load data and check channel identity
    cFile = [opts.fPath filesep 'Analog_' num2str(trialNr) '.dat']; %current file to be read
    [~,Analog] = Widefield_LoadData(cFile,'Analog'); %load analog data
    Analog = double(Analog);
    
    offset = 0;
    bFrame = find(temp < min(temp)*.75); %index for black frames
    if bFrame(1) == 1 %if first frame is dark, remove initial frames from data until LEDs are on
        %remove initial dark frames
        offset = find(diff(bFrame) > 1, 1);
        data(:,:,1:offset) = [];
        dSize = size(data);
        temp(1:offset) = [];
        frameTimes(1:offset) = [];
    end
    
    % check if pre- and poststim are given. use all frames if not.
    if ~isfield(opts,'preStim') || ~isfield(opts,'preStim')
        opts.preStim = 0;
        opts.postStim = inf;
    else
        opts.preStim = ceil(opts.preStim * sRate);
        opts.postStim = ceil(opts.postStim * sRate);
    end
    
    if any(~isnan(opts.trigLine)) || any(opts.trigLine > size(Analog,1))
        
        trace = Analog(opts.trigLine,:); %blue and violet light trigger channels
        trace = zscore(trace(1,end:-1:1) - trace(2,end:-1:1)); %invert and subtract to check color of last frame
        trace(round(diff(trace)) ~= 0) = 0; %don't use triggers that are only 1ms long
        lastBlue = find(trace > 1, 1);
        lastHemo = find(trace < -1,1);
        
        blueLast = lastBlue < lastHemo;
        if isempty(lastBlue) || isempty(lastHemo)
            warning(['Failed to find trigger signals. lastBlue: ' num2str(lastBlue) '; lastHemo: ' num2str(lastHemo) '; trialNr: ' num2str(trialNr)])
        end
        
        bFrame = find(temp < min(temp)*.75); %index for first black frame (assuming the last frame is really a dark frame)
        bFrame(bFrame < round(size(temp,2) / 2)) = []; %make sure black frame is in the second half of recording.
        bFrame = bFrame(1);
        blueFrames = false(1,length(temp));
        
        if blueLast %last frame before black is blue
            if rem(bFrame,2) == 0 %blue frames (bFrame - 1) have uneven numbers
                blueFrames(1:2:dSize(end)) = true;
            else
                blueFrames(2:2:dSize(end)) = true;
            end
            lastFrame = size(Analog,2) - lastBlue; %index for end of last frame
        else %last frame before black is violet
            if rem(bFrame,2) == 0 %blue frames (bFrame - 2) have even numbers
                blueFrames(2:2:dSize(end)) = true;
            else
                blueFrames(1:2:dSize(end)) = true;
            end
            lastFrame = size(Analog,2) - lastHemo; %index for end of last frame
        end
        
        frameTimes = (frameTimes - frameTimes(bFrame - 1)) + lastFrame; %realign frameTime based on time of last non-dark frame
        blueFrames = blueFrames(frameTimes < size(Analog,2));
        blueFrames(bFrame - 1:end) = []; %exclude black and last non-black frame
        
        blueTimes = frameTimes(blueFrames);
        hemoTimes = frameTimes(~blueFrames);
        
        blueData = data(:,:,blueFrames);
        hemoData = data(:,:,~blueFrames);
        
        % find all triggers in stim line and choose the one that is on the longest as the true stimulus trigger
        % apparently this line can be contaminated by noise
        stimOn = find(diff(double(Analog(opts.stimLine,:)) > 1500) == 1);
        ind = find((find(diff(double(Analog(opts.stimLine,:)) > 1500) == -1) - stimOn) > 2,1); %only use triggers that are more than 2ms long
        stimOn = stimOn(ind) + 1;
        
        if ~isinf(opts.postStim)
            if isempty(stimOn) || isempty(find((blueTimes - stimOn) > 0, 1))
                error(['No stimulus trigger found. Current file: ' cFile])
            else
                blueStim = find((blueTimes - stimOn) > 0, 1);  %first blue frame after stimOn
                if blueStim <= opts.preStim + 1 %if stim occured earlier as defined by preStim
                    blueStim = opts.preStim + 1 ; %use earliest possible frame for stimulus onset
                    if (find(hemoTimes - blueTimes(opts.preStim + 1) > 0, 1) - 1) <= opts.preStim %make sure there are enough frames for hemo channel
                        blueStim = blueStim + 1 ; %use earliest possible frame for stimulus onset
                    end
                    fprintf('Warning: StimOn is too early. Starting at %d instead of %d to provide enough baseline frames.\n',blueStim, find((blueTimes - stimOn) > 0, 1))
                    falseAlign = true;
                end
            end
            
            if length(blueTimes) < (blueStim + opts.postStim - 1)
                lastBlueIdx = length(blueTimes);
            else
                lastBlueIdx = blueStim + opts.postStim - 1;
            end
            
            hemoStim = find(hemoTimes - blueTimes(blueStim) > 0, 1) - 1;  %first hemo frame before blueFrame after stimOn :)
            if length(hemoTimes) < (hemoStim + opts.postStim - 1)
                lastHemoIdx = length(hemoTimes);
            else
                lastHemoIdx = hemoStim + opts.postStim - 1;
            end
            
            %make sure both channels have equal length
            chanDiff = length(blueStim - opts.preStim : lastBlueIdx) - length(hemoStim - opts.preStim : lastHemoIdx);
            if chanDiff < 0 %less blue frames, cut hemo frame
                lastHemoIdx = lastHemoIdx + chanDiff;
            elseif chanDiff > 0 %less hemo frames, cut blue frame
                lastBlueIdx = lastBlueIdx - chanDiff;
            end
            
            %check if enough data is remaining after stimulus onset
            if lastBlueIdx - blueStim + 1  < opts.postStim
                fprintf('Warning: StimOn is too late. Using %d instead of %d frames for poststim data.\n',lastBlueIdx - blueStim, opts.postStim)
                falseAlign = true;
            end
            
            blueTimes = blueTimes(blueStim - opts.preStim : lastBlueIdx); %get blue frame times before and after stimOn
            blueData = blueData(:, :, blueStim - opts.preStim : lastBlueIdx); %get blue frame data before and after stim on
            
            hemoTimes = hemoTimes(hemoStim - opts.preStim : lastHemoIdx); %get hemo frame times before and after stimOn
            hemoData = hemoData(:,:,hemoStim - opts.preStim : lastHemoIdx); %get hemo frame data before and after stim on
        else
            %make sure both channels have equal length
            chanDiff = size(blueData,3) - size(hemoData,3);
            if chanDiff < 0 %less blue frames, cut hemo frame
                hemoData = hemoData(:,:,1:end+chanDiff);
                hemoTimes = hemoTimes(1:end+chanDiff);
            elseif chanDiff > 0 %less hemo frames, cut blue frame
                blueData = blueData(:,:,1:end-chanDiff);
                blueTimes = blueTimes(1:end-chanDiff);
            end
            
            % check for stimulus onset
            if isempty(stimOn) || isempty(find((blueTimes - stimOn) > 0, 1))
                warning(['No stimulus trigger found. Current file: ' cFile])
                stimOn = NaN;
            else
                stimOn = find((blueTimes - stimOn) > 0, 1);  %first blue frame after stimOn
            end
        end
        
        %make sure each frame has a corresponding trigger signal
        frameDiff = floor(median(diff(frameTimes)));
        frameJitter = sum(diff(frameTimes) > (frameDiff + 1));
        frameDiff = frameDiff - (1 - rem(frameDiff,2)); %make sure, inter-frame interval is uneven number
        blueInd = bsxfun(@minus, repmat(round(blueTimes),1,frameDiff),-floor(frameDiff/2):floor(frameDiff/2))'; %get index for each frameTime and get half the IFI before and after frame was acquired to check trigger signal.
        hemoInd = bsxfun(@minus, repmat(round(hemoTimes),1,frameDiff),-floor(frameDiff/2):floor(frameDiff/2))'; %get index for each frameTime and get half the IFI before and after frame was acquired to check trigger signal.
        
        blueTrig = zscore(Analog(opts.trigLine(1), blueInd(:))); %blue light trigger channel
        blueTrig = reshape(blueTrig, [], length(blueTimes));
        blueTrig = any(zscore(blueTrig) > 0.5);
        
        hemoTrig = zscore(Analog(opts.trigLine(2), hemoInd(:))); %blue light trigger channel
        hemoTrig = reshape(hemoTrig, [], length(hemoTimes));
        hemoTrig = any(zscore(hemoTrig) > 0.5);
        
        if sum(blueTrig) ~= length(blueTimes) || sum(hemoTrig) ~= length(hemoTimes)
            disp(['Potential frame time violations: ' num2str(frameJitter)])
            disp(['Confirmed blue trigs: ' num2str(sum(blueTrig)) '; blue frame times: ' num2str(length(blueTimes))])
            disp(['Confirmed hemo trigs: ' num2str(sum(hemoTrig)) '; hemo frame times: ' num2str(length(hemoTimes))])
            
            if frameJitter < (length(blueTimes) - sum(blueTrig)) || frameJitter < (length(hemoTimes) - sum(hemoTrig)) %more unaccounted triggers as might be explained by frame time jitter
                falseAlign = true; %potential for misaligned channel separation. Don't use trial.
                disp('Flagged file for rejection.')
            end
            disp(['Current file: ' cFile])
        end
    else
        blueData = data(:,:,1 : opts.preStim + opts.postStim);
        blueTimes = frameTimes(1: opts.preStim + opts.postStim);
        hemoData = NaN;
        hemoTimes = NaN;
        falseAlign = true;
        stimOn = [];
    end
else
    data(:,:,1:offset) = [];
    blueData = data(:,:,blueFrames);
    hemoData = data(:,:,~blueFrames);
    blueTimes = [];
    hemoTimes = [];
    stimOn = [];
end
end


function [header,data] = Widefield_LoadData(cPath,Condition,dataType,pInd)
% short routine to load data from WidefieldImager code.
% cPath is the path of the file that should be opened. Condition is the
% type of data file which will determine the way the data file is read.
% Optional input 'pInd' defines a single pixel from which a data vector should
% be extracted. pInd is a two-value vector for x-y pixel coordinates. This
% means X and Y for image coordinates NOT matlab coordinates (which are
% usually inverse).

if ~exist('dataType','var') || isempty(dataType)
    dataType = 'uint16';
end

if ~exist('pInd','var') || isnan(pInd)
    pInd = [];
else
    if length(pInd) ~= 2 || ~isnumeric(pInd)
        error('Invalid input for index of selected pixel')
    end
end

fID = fopen(cPath);
switch lower(Condition)
    case 'analog'
        hSize = fread(fID,1,'double'); %header size
        header = fread(fID,hSize,'double'); %Metadata. Default is: 1 = Time of Acquisition onset, 2 = Number of channels, 3 = number of values per channel
        data = fread(fID,[header(end-1),header(end)],[dataType '=>' dataType]); %get data. Last 2 header values should contain the size of the data array.
    case 'frames'
        hSize = fread(fID,1,'double'); %header size
        header = fread(fID,hSize,'double'); %Metadata. Default is: 1:x = Absolute timestamps for each frame, Last 4 values: Size of each dimensions in the matrix
        
        if ~isempty(pInd) %if extracting single pixel information
            imSize = (header(find(diff(header) < -1e3) + 1)*header(find(diff(header) < -1e3) + 2))-1; %number of datapoints to make a single image minus one. skip that many datapoints to stick to the same pixel when using fread.
            imStart = ((pInd(1)-1)*header(find(diff(header) < -1e3) + 1))+pInd(2)-1; %first value for selected pixel
            fseek(fID,imStart*2,'cof'); %shift file pointer to the right pixel to start data extraction from file
            data = fread(fID,header(find(diff(header) < -1e3) + 4),[dataType '=>' dataType],imSize*2); %get data.
            if length(data) ~= header(end)
                error('Could not extract all data values from pixel')
            end
        else
            data = fread(fID,[prod(header(find(diff(header) < -1e3) + 1 : end)),1],[dataType '=>' dataType]); %get data. Last 4 header values should contain the size of the data array.
            if length(data) ~= prod(header(find(diff(header) < -1e3) + 1 : end)) %if insufficient data is found in .dat file. Sometimes fread does not get all values from file when reading from server.
                fclose(fID);fID = fopen(cPath); %try loading data again
                hSize = fread(fID,1,'double'); %header size
                header = fread(fID,hSize,'double'); %Metadata. Defautlt is: 1:x = Absolute timestamps for each frame, Last 4 values: Size of each dimensions in the matrix
                data = fread(fID,[prod(header(find(diff(header) < -1e3) + 1 : end)),1],[dataType '=>' dataType]); %get data. Last 4 header values should contain the size of the data array.
            end
            data = reshape(data,header(find(diff(header) < -1e3) + 1 : end)'); %reshape data into matrix
        end
end
fclose(fID);
end
function splitChannelsTemp(opts,trialNr)
% Code to separate blue and violet channel from widefield data. This needs
% analog data that contains a stimulus onset from which a pre and
% poststimulus dataset can be returned. Also requires trigger channels for
% blue/violet LEDs and some dark frames at the end.
% This code is the newer version of Widefield_CheckChannels.

falseAlign = false;

%% load data and check channel identity
cFile = [opts.fPath filesep 'Analog_' num2str(trialNr) '.dat']; %current file to be read
[~,Analog] = loadRawData(cFile,'Analog'); %load analog data
Analog = double(Analog);

load([opts.fPath filesep 'frameTimes_' num2str(trialNr, '%04i') '.mat'], 'imgSize', 'frameTimes'); %get data size
frameTimes = frameTimes * 86400*1e3; %convert to seconds
cFile = [opts.fPath filesep opts.fName '_' num2str(trialNr, '%04i') '.dat']; %current file to be read
[~, data] = loadRawData(cFile,'Frames',[], imgSize); %load video data

%reshape data to compute mean frame intensities
dSize = size(data);
data = reshape(data,[],dSize(end));
temp = zscore(mean(single(data)));
data = squeeze(reshape(data,dSize));

bFrame = find(temp < min(temp)*.75); %index for black frames
if bFrame(1) == 1 %if first frame is dark, remove initial frames from data until LEDs are on
    %remove initial dark frames
    cIdx = find(diff(bFrame) > 1, 1); 
    data(:,:,1:cIdx) = [];
    dSize = size(data);
    temp(1:cIdx) = [];
    frameTimes(1:cIdx) = [];
end

%determine imaging rate - either given as input or determined from data
if isfield(opts,'frameRate')
    sRate = opts.frameRate;
else
    sRate = 1000/(mean(diff(frameTimes))*2);
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
   
    bFrame = find(temp < min(temp)*.75); %index for first black frame (assuming the last frame is really a dark frame)
    bFrame(bFrame < round(size(temp,2) / 2)) = []; %make sure black frame is in the second half of recording.
    
    if isempty(bFrame); bFrame = length(temp); end %if there are no black frames, just use the last one
    bFrame = bFrame(1);
    
    %get number of rejected frames before data was saved
    nrFrames = bFrame - 1; %number of exposed frames in the data
    nrTriggers = sum(diff(smooth(Analog(8,:),10) > 1000) == 1) + sum(diff(smooth(Analog(9,:),10) > 1000) == 1); %nr of triggers
    removedFrames = nrTriggers - nrFrames;
    save([opts.fPath filesep 'frameTimes_' num2str(trialNr, '%04i') '.mat'], 'imgSize', 'frameTimes', 'removedFrames'); %get data size

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
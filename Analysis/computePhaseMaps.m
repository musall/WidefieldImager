function computePhaseMaps(fPath,Animal,nTrials,winSize,smth,rotateImage)
dbstop Widefield_ComputePhaseMaps 113;

%% Set basic variables
fName = 'Vc'; %name of data file

if ~strcmpi(fPath(end),filesep)
    fPath = [fPath filesep];
end

if ~exist('nTrials','var')
    nTrials = 1;                        %number of trials for phasemaps. If nCycles <=0 or > the number of cyles in the data. It will use all available cycles for a single phasemap.
end

if ~exist('winSize','var')
    winSize = 10;                         %size of the field in mm. This is to determine the spatial binning to get as close to 40pix/mm as possible. This is advised for the visual segmentation code.
end

if ~exist('smth','var')
    smth = 1.5;                         %smoothing factor for gaussian smooth on maps
end

if ~exist('rotateImage','var')
    rotateImage = 0;                    %rotate image if the orientation is not as required
end


%% data source
disp(['Current path: ' fPath]); tic

%Load svd-compressed data. Contains all imaging data and timestamps in ms for each frame.
load([fPath fName], 'Vc', 'U', 'frameCnt', 'stimTime');

%load settings
cFile = ls([fPath Animal '*settings.mat']);
load([fPath cFile], 'StimData');
Trials = 1 : str2double(StimData.handles.NrTrials);
sRate = round(1 / StimData.sRate) / 2; %frame rate in Hz

%% get different bar direction and orentiations for individual trials
BarOrient = StimData.VarVals(strcmpi(StimData.VarNames,'BarOrient'),:); %get bar orientation for all trials
BarDirection = StimData.VarVals(strcmpi(StimData.VarNames,'BarDirection'),:); %get bar direction for all trials

iBarConds(1,:) = BarOrient == 1  & BarDirection == 0; % horizontal moving down
iBarConds(2,:) = BarOrient == 1  & BarDirection == 1; % horizontal moving up
iBarConds(3,:) = BarOrient == 0  & BarDirection == 0; % vertial, moving left to right
iBarConds(4,:) = BarOrient == 0  & BarDirection == 1; % vertical, moving right to left

% get duration and bar speed for individual trials. this is assumed to be constant by now.
StimDur = StimData.VarVals(strcmpi(StimData.VarNames,'StimDuration') | strcmpi(StimData.VarNames,'trialDuration'),:); %Duration of a given trial
barFreq = StimData.VarVals(strcmpi(StimData.VarNames,'cyclesPerSecond'),:); %bar speed in a given trial
numCycles = unique(StimDur(Trials)./(1./barFreq(Trials))); %number of cycles in current trial

% check for number of cycles and trials in dataset
if length(numCycles) ~= 1
    error('Number of cycles per trial is inconsistent. This code is not meant to handle that.')
end
disp(['Trials per condition: ' num2str(length(Trials)/4) ' - Computing phasemaps from [' num2str(nTrials) '] trials']);
TrialCnts = ones(4,length(nTrials));
avgData = cell(4,length(nTrials)); % averaged sequence for each condition and trialcount
fTransform = cell(4,length(nTrials)); % fourier transforms for each condition and trialcount 

%% get single trials and average for each condition to end up with one sequence
condCnt = zeros(1,4); %counter for how many sequences were saved in each condition
for iTrials = Trials
    %% load data
    trialStart = sum(frameCnt(2,1:iTrials)) - frameCnt(2,iTrials) + 1; %first frame of current trial
    Data = reshape(U,[],size(U,3)) * Vc(:,trialStart : sum(frameCnt(2,1:iTrials)));
    Data = reshape(Data, size(U,1), size(U,2), []);  
    
    %% compute the amount of required frames and collect from data
    stimFrames = StimDur(iTrials)*sRate; %required frames for stim sequence
    Data = Data(:,:,stimTime(iTrials): stimTime(iTrials) + stimFrames - 1); %only use stimulation part of the data;

    binSize = floor((max(size(Data(:,:,1)))/winSize)/40); %compute binsize to get closest to 40 pixels/mm (recommended for segmentation code).
    if winSize > 1 && winSize < inf
        bData = arrayResize(Data,binSize); %do spatial binning
    else
        bData = Data;
    end
    clear Data

    %% running average
    for x = 1:length(nTrials)
        cTrialCnt = rem(condCnt(iBarConds(:,iTrials))+1,nTrials(x)); %current trial cycle for running average (reset to 1, when 'nTrials' is reached)
        
        if cTrialCnt == 1 || condCnt(iBarConds(:,iTrials)) == 0 %start cycle for running average
            avgData{iBarConds(:,iTrials),x} = bData; %starting dataset for running average with set trialcount
        else
            avgData{iBarConds(:,iTrials),x} = (avgData{iBarConds(:,iTrials),x}.*cTrialCnt + bData) ./ (cTrialCnt+1); %produce running average
        end
        
        if cTrialCnt == 0 %reached requested trialcount. Compute fourier transform and increase counter
            temp = fft(avgData{iBarConds(:,iTrials),x},[],3);
            fTransform{iBarConds(:,iTrials),x}(TrialCnts(iBarConds(:,iTrials),x),:,:) = squeeze(temp(:,:,numCycles+1)); clear temp
            TrialCnts(iBarConds(:,iTrials),x) = TrialCnts(iBarConds(:,iTrials),x)+1; %increase counter for running average when required trialcount is reached
        end
    end
    condCnt(iBarConds(:,iTrials)) =  condCnt(iBarConds(:,iTrials))+1;
    if rem(iTrials, 10) == 0
        disp(['Done loading trial ' int2str(iTrials) '/' int2str(max(Trials))]);
    end
end
clear Data bData

%% do fft analysis to get phase and magnitude maps
for iTrials = 1
    Cnt = 1;
    for iConds = [1 3] %this expects 4 directions to construct horizontal and vertical map
        for iRuns = 1:size(fTransform{iConds,iTrials},1)
        
           magMaps{Cnt,iRuns} = imrotate(squeeze(abs(fTransform{iConds,iTrials}(iRuns,:,:).*fTransform{iConds+1,iTrials}(iRuns,:,:))),rotateImage); %combined magnitude map.
           phaseMaps{Cnt,iRuns} = imrotate(squeeze(angle(fTransform{iConds,iTrials}(iRuns,:,:)./fTransform{iConds+1,iTrials}(iRuns,:,:))),rotateImage); %combined phase map (horizontal is iConds=1, vertical is iConds =3)
           phaseMaps{Cnt,iRuns}(isnan(phaseMaps{Cnt,iRuns}(:))) = 0;

        end
        cPhaseMaps{Cnt,iTrials} = median(cat(3,phaseMaps{Cnt,:}),3);
        Cnt = Cnt+1;
    end
    cMagMaps{iTrials} = median(cat(3,magMaps{:}),3);
    cMagMaps{iTrials} =(cMagMaps{iTrials}-min(cMagMaps{iTrials}(:)))./(max(cMagMaps{iTrials}(:))- min(cMagMaps{iTrials}(:))); %normalize between 0 and 1

    % compute visual field sign maps. First compute gradients and atan - same as in the 2014 Callaway paper.
    for iRuns = 1:size(fTransform{iConds,iTrials},1)
        [dhdx, dhdy] = gradient(phaseMaps{1,iRuns});
        [dvdx, dvdy] = gradient(phaseMaps{2,iRuns});
        
        graddir_hor = atan2(dhdy,dhdx);
        graddir_vert = atan2(dvdy,dvdx);
        vdiff = exp(1i*graddir_hor) .* exp(-1i*graddir_vert);
        
        VFS{iRuns} = sin(angle(vdiff)); %Visual field sign map
        VFS{iRuns} = spatialFilterGaussian(VFS{iRuns},smth);
    end
    
    cVFS{1,iTrials} = median(cat(3,VFS{:}),3);
    clear magMaps phaseMaps VFS
    
    h = figure;
    subplot(2,2,1);
    imagesc(cPhaseMaps{1,iTrials});axis image; colormap hsv; colorbar; freezeColors;
    title(['Horizontal - nTrials = ' num2str(nTrials(iTrials))]);
    subplot(2,2,2);
    imagesc(cPhaseMaps{2,iTrials});axis image; colormap hsv; colorbar; freezeColors;
    title(['Vertical - nTrials = ' num2str(nTrials(iTrials))]);
    subplot(2,2,3);
    imagesc(cMagMaps{iTrials});axis image; colorbar; colormap jet;
    title('Mean Magnitude');
    subplot(2,2,4);
    imagesc(spatialFilterGaussian(cVFS{1,iTrials},smth)); axis image;colorbar
    caxis([-0.5 0.5])
    title(['VisualFieldSign - binSize = ' num2str(binSize) '; smth = ' num2str(smth)]);
    savefig(h,[fPath '\' Animal '_phaseMap_allPlots_ ' int2str(nTrials(iTrials)) '_trials.fig']);
    h.PaperUnits = 'inches';
    set(h, 'PaperPosition', [0 0 15 15]);
    saveas(h,[fPath '\' Animal '_phaseMap_allPlots_ ' int2str(nTrials(iTrials)) '_trials.jpg'])
    clear h
end

trialSelect = 1; %use this to select which trialcount should be used for subsequent figures


%% get phase and amplitude + vessel map for plotting
plotPhaseMap = spatialFilterGaussian(cVFS{1,trialSelect},smth);
plotPhaseMap = imresize(plotPhaseMap,binSize);

plotAmpMap = spatialFilterGaussian(imresize(cMagMaps{trialSelect},binSize),25); %smoothed magnitude map
plotAmpMap =(plotAmpMap-min(plotAmpMap(:)))./(max(plotAmpMap(:))- min(plotAmpMap(:))); %normalize between 0 and 1
plotAmpMap = spatialFilterGaussian(plotAmpMap,25); %smoothed magnitude map

%get vessel image.
cFile = ls([fPath 'Snapshot_1.mat']);
load([fPath cFile],'snap');
snap = double(imrotate(snap,rotateImage));
snap =(snap-min(snap(:)))./(max(snap(:))- min(snap(:))); %normalize between 0 and 1

if size(snap,1) > size(plotPhaseMap,1) %resize if vessel map is larger
    plotPhaseMap = imresize(plotPhaseMap,size(snap,1)/size(plotPhaseMap,1));
    plotAmpMap = imresize(plotAmpMap,size(snap,1)/size(plotAmpMap,1));
end

%% plot vessel map with overlayed sign map
h = figure;
imagesc(plotPhaseMap); colormap jet; 
caxis([-0.5 0.5])
title([Animal ' - PhaseMap']);axis image
savefig(h,[fPath Animal '_RawPhaseMap.fig']);
saveas(h,[fPath Animal '_RawPhaseMap.jpg']);
save([fPath 'plotPhaseMap.mat'],'plotPhaseMap');
save([fPath 'plotAmpMap.mat'],'plotAmpMap');
save([fPath 'cMagMaps.mat'],'cMagMaps');
save([fPath 'cPhaseMaps.mat'],'cPhaseMaps');
save([fPath 'cVFS.mat'],'cVFS');

h = figure;
imagesc(snap);axis image; colormap gray; freezeColors; hold on
vfsIm = imagesc(plotPhaseMap); colormap jet; 
caxis([-0.5 0.5]);
set(vfsIm,'AlphaData',plotAmpMap); axis image
title([Animal ' - PhaseMap + Vesselmap'])
savefig(h,[fPath Animal '_phaseMap.fig']);
saveas(h,[fPath Animal '_phaseMap.jpg'])



%nested functions
function img = spatialFilterGaussian(img, sigma)
if sigma > 0 && (numel(img) ~=  sum(sum(isnan(img))))
    hh = fspecial('gaussian',size(img),sigma);
    hh = hh/sum(hh(:));
    img = ifft2(fft2(img).*abs(fft2(hh)));
end

function [header,data] = Widefield_LoadData(cPath,Condition,varargin)
% short routine to load data from WidefieldImager code.
% cPath is the path of the file that should be opened. Condition is the
% type of data file which will determine the way the data file is read.
% Optional input 'pInd' defines a single pixel from which a data vector should
% be extracted. pInd is a two-value vector for x-y pixel coordinates. This
% means X and Y for image coordinates NOT matlab coordinates (which are
% usually inverse).

if ~isempty(varargin)
    pInd = varargin{1}; %index for pixel to be extracted
    if length(pInd) ~= 2 || ~isnumeric(pInd)
        error('Invalid input for index of selected pixel')
    end
end

fID = fopen(cPath);
switch lower(Condition)
    case 'analog'
        hSize = fread(fID,1,'double'); %header size
        header = fread(fID,hSize,'double'); %Metadata. Default is: 1 = Time of Acquisition onset, 2 = Number of channels, 3 = number of values per channel
        data = fread(fID,[header(end-1),header(end)],'uint16=>uint16'); %get data. Last 2 header values should contain the size of the data array.
    case 'frames'
        hSize = fread(fID,1,'double'); %header size
        header = fread(fID,hSize,'double'); %Metadata. Default is: 1:x = Absolute timestamps for each frame, Last 4 values: Size of each dimensions in the matrix
        
        if ~isempty(varargin) %if extracting single pixel information
            imSize = (header(find(diff(header) < -1e4) + 1)*header(find(diff(header) < -1e4) + 2))-1; %number of datapoints to make a single image minus one. skip that many datapoints to stick to the same pixel when using fread.
            imStart = ((pInd(1)-1)*header(find(diff(header) < -1e4) + 1))+pInd(2)-1; %first value for selected pixel
            fseek(fID,imStart*2,'cof'); %shift file pointer to the right pixel to start data extraction from file
            data = fread(fID,header(find(diff(header) < -1e4) + 4),'uint16=>uint16',imSize*2); %get data.
            if length(data) ~= header(end)
                error('Could not extract all data values from pixel')
            end
        else
            data = fread(fID,[prod(header(find(diff(header) < -1e4) + 1 : end)),1],'uint16=>uint16'); %get data. Last 4 header values should contain the size of the data array.
            if length(data) ~= prod(header(find(diff(header) < -1e4) + 1 : end)) %if insufficient data is found in .dat file. Sometimes fread does not get all values from file when reading from server.
                fclose(fID);fID = fopen(cPath); %try loading data again
                hSize = fread(fID,1,'double'); %header size
                header = fread(fID,hSize,'double'); %Metadata. Defautlt is: 1:x = Absolute timestamps for each frame, Last 4 values: Size of each dimensions in the matrix
                data = fread(fID,[prod(header(find(diff(header) < -1e4) + 1 : end)),1],'uint16=>uint16'); %get data. Last 4 header values should contain the size of the data array.
            end
            data = reshape(data,header(find(diff(header) < -1e4) + 1 : end)'); %reshape data into matrix
        end
end
fclose(fID);
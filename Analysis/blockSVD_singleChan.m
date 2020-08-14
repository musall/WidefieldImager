function [bV, bU, blockInd, frameCnt, blueAvg] = blockSVD_singleChan(opts)
% Code to convert raw data from a given imaging experiment to a
% low-dimensional representation. Imaging data consists of a single imaging
% channel and all data is used.
% opts.fPath denotes the path to the example data folder, containing raw video 
% data as .dat files and according timestamps for each frame. 
% nrBlocks is the total number of blocks used for the svd. Sqrt of blocks
% has to be an even number - nrBlocks is rounded down if needed.
% overlap determines the number of pixels with which individual blocks are
% overlapping to avoid edge effects.
%
% Inputs: 
% opts: A structure that contains option for dimensionality reduction. 
%       Refer to 'Tutorial_dimReduction' for an example.
%
% Outputs:
% bV: A cell array that contains temporal components for each block over
%     all frames in the session.
% bU: A cell array that contains spatial components for each block.
%     Convolve with bV to reconstruct original imaging data.
% blockInd: Index that can be used to place spatial components in 'bU' in a
%           larger larger matrix that is the size of the original widefield 
%           data. Use 'blueAvg' or 'hemoAvg' to determine the size of a
%           frame in the widefield data.
% frameCnt: Indicates the number of frames in each trial. First ror is the
%           trial number, second row is number of frames.
% blueAvg: Mean image of the imaging data. Can be used to get vessel image
%          or determine size of a single frame in the widefield data.
%
% Simon Musall, 8/11/2020

%% check options
if ~strcmpi(opts.fPath(end),filesep)
    opts.fPath = [opts.fPath filesep];
end

% these should be provided as inputs
nrBlocks = opts.nrBlocks;
overlap = opts.overlap;
opts.plotChans = false; %plot results for blue and violet channel when loading raw data. Will make the code break on the HPC.
opts.verbosity = false; %flag to supress warnings from 'splitChannel' code.

if nrBlocks ~= floor(sqrt(nrBlocks))^2
    fprintf('Chosen nrBlocks (%d) cannot be squared. Using %d instead.\n', nrBlocks, floor(sqrt(nrBlocks))^2)
    nrBlocks = floor(sqrt(nrBlocks))^2;
end

disp('===============');
disp(opts.fPath); 
fprintf('Sampling rate: %dHz; Using %d blocks for SVD\n', opts.frameRate,nrBlocks);
disp(datestr(now));
tic;

%% check if rawData files are present and get trialnumbers for all files
rawCheck = dir([opts.fPath filesep opts.fName '*dat']);

fileCnt = size(rawCheck,1);
for iFiles = 1 : fileCnt %go through files and check for trialNr. this is used for transformation to .mat files later.
    temp = textscan(rawCheck(iFiles).name,'%s%f%s','Delimiter','_');
    trials(iFiles) = temp{2};
end
trials = sort(trials);

%% get reference images for motion correction
cFile = [opts.fPath filesep opts.fName '_' num2str(trials(1)) '.dat']; %first file to be read
[~, blueData] = Widefield_LoadData(cFile,'Frames'); %load video data
        
blueData = single(squeeze(blueData));
blueRef = fft2(median(blueData,3)); %blue reference for alignment
save([opts.fPath 'blueRef.mat'],'blueRef');

if opts.useGPU
    blueRef = gpuArray(blueRef);
end

%% get index for individual blocks
indImg = reshape(1:numel(blueRef),size(blueRef)); %this is an 'image' with the corresponding indices
blockSize = ceil((size(blueRef) + repmat(sqrt(nrBlocks) * overlap, 1, 2))/sqrt(nrBlocks)); %size of each block
blockInd = cell(1, nrBlocks);

Cnt = 0;
colSteps = (0 : blockSize(1) - overlap : size(blueRef,1)) + 1; %steps for columns
rowSteps = (0 : blockSize(2) - overlap : size(blueRef,2)) + 1; %steps for rows
for iRows = 1 : sqrt(nrBlocks)
    for iCols = 1 : sqrt(nrBlocks)
        
        Cnt = Cnt + 1;
        % get current block and save index as vector
        colInd = colSteps(iCols) : colSteps(iCols) + blockSize(1) - 1; 
        rowInd = rowSteps(iRows) : rowSteps(iRows) + blockSize(2) - 1;
        
        colInd(colInd > size(blueRef,1)) = [];
        rowInd(rowInd > size(blueRef,2)) = [];
        
        cBlock = indImg(colInd, rowInd);
        blockInd{Cnt} = cBlock(:);
        
    end
end
save([opts.fPath 'blockInd.mat'],'blockInd');

%% perform image alignement for separate channels and collect data in mov matrix
blueAvg = zeros([size(blueData,1), size(blueData,2), fileCnt],'uint16'); %average for mean correction. Collect single session averages to verify correct channel separation.
blueFrameTimes = cell(1, fileCnt);
if ~exist([opts.fPath 'blockData'], 'dir')
    mkdir([opts.fPath 'blockData']);
end

frameCnt = NaN(2,fileCnt, 'single'); %use thise to report how many frames were collected in each trial
blueBlocks = zeros(1,nrBlocks);
for iTrials = 1:fileCnt
    
    cFile = [opts.fPath filesep opts.fName '_' num2str(trials(iTrials)) '.dat']; %first file to be read
    [header, blueData] = Widefield_LoadData(cFile,'Frames'); %load video data
    blueTimes = header(1:size(blueData,3)); %get frametimes from header

    frameCnt(1,iTrials) = trials(iTrials); %trialNr
    frameCnt(2,iTrials) = size(blueData,3); %nr of frames
    
    if opts.useGPU
        blueData = gpuArray(blueData);
    end
    
    %perform image alignment for both channels
    for iFrames = 1:size(blueData,3)
        [~, temp] = dftregistration(blueRef, fft2(blueData(:, :, iFrames)), 10);
        blueData(:, :, iFrames) = abs(ifft2(temp));
    end
    blueData = gather(blueData);
    
    % keep baseline average for each trial
    blueAvg(:,:,iTrials) = mean(blueData(:,:,1:opts.baselineFrames),3);
    
    %keep timestamps for all frames
    blueFrameTimes{iTrials} = blueTimes;

    if rem(iTrials,10) == 0
        fprintf(1, 'Loading session %d out of %d\n', iTrials,length(trials));
    end
    
    % save data in individual blocks. single file for each trial/block. Will delete those later.
    blueData = reshape(blueData, [], size(blueData,3));
    for iBlocks = 1:nrBlocks
        if iTrials == 1
            blueBlocks(iBlocks) = fopen([opts.fPath 'blockData' filesep 'blueBlock' num2str(iBlocks) '.dat'], 'Wb');
        end

        bBlock = blueData(blockInd{iBlocks}, :);
        fwrite(blueBlocks(iBlocks), bBlock,'uint16'); %write data from current trial to block file
        
        if iTrials == fileCnt
            fclose(blueBlocks(iBlocks));
        end
    end
    clear blueData blueTimes
end
disp('Binary files created!'); toc;

clear blueRef 
save([opts.fPath 'trials.mat'],'trials'); %save trials so order of analysis is consistent

%save frametimes for blue/hemo trials
save([opts.fPath 'blueFrameTimes.mat'],'blueFrameTimes', 'trials');

%save averages in case you need them later
save([opts.fPath 'blueAvg.mat'],'blueAvg');

%take average over all trials for subsequent mean correction
blueAvg = mean(single(blueAvg),3);

%% compress each block with SVD
bU = cell(nrBlocks,1); bV = cell(nrBlocks,1);
for iBlocks = 1 : nrBlocks
    
    % load current block
    fID = fopen([opts.fPath 'blockData' filesep 'blueBlock' num2str(iBlocks) '.dat'], 'r');
        
    allBlock = fread(fID, 'uint16'); fclose(fID(1));
    allBlock = reshape(allBlock, size(blockInd{iBlocks},1), size(cat(1,blueFrameTimes{:}),1))';
    delete([opts.fPath 'blockData' filesep 'blueBlock' num2str(iBlocks) '.dat']);

    % run SVD on current block
    allBlock = reshape(allBlock,size(allBlock,1),[])'; %combine channels and transpose (this is faster if there are more frames as pixels)
    [bV{iBlocks}, s, bU{iBlocks}] = fsvd(allBlock,opts.blockDims); %U and V are flipped here because we transpoed the input.
    bV{iBlocks} = gather(s * bV{iBlocks}'); %multiply S into V, so only U and V from here on
    bU{iBlocks} = gather(bU{iBlocks});
    clear allBlock
    
    if rem(iBlocks, round(nrBlocks / 5)) == 0
        fprintf(1, 'Converting block %d out of %d\n', iBlocks, nrBlocks);
    end
end

% save blockwise SVD data from both channels
save([opts.fPath 'bV.mat'], 'bU', 'bV', 'blockInd', 'opts', '-v7.3');
disp('Blockwise SVD complete'); toc;

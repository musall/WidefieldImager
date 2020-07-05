%construct path to data folder and give some basic info
opts.fPath = [pwd filesep 'DemoRec' filesep]; %path to demo recording
opts.fName = 'Frames'; %name of imaging data files.
opts.nrBlocks = 49; %nr of blocks for svd
opts.overlap = 20; % pixel overlap between blocks
opts.dimCnt = 200; %nr of components in the final dataset
opts.blockDims = 25; %number of dimensions from SVD per block
opts.stimLine = 3; %analog line that contains stimulus trigger.
opts.trigLine = [6 7]; %analog lines for blue and violet light triggers.
opts.useGPU = false; %flag to use GPU acceleration

% check for handles file from WidefieldImager to get frame rate.
try
    load([opts.fPath 'handles.mat'],'handles'); %get handles from widefieldimager
    opts.frameRate = str2double(handles.FrameRate.String)/2; %frame rate of individual channels. With alternating illumination this is the absolute frameRate / 2.
    close(handles.WidefieldImager); %loading handles bring sup the GUI so close that again.
catch
    opts.frameRate = 15; %single-channel sampling rate in Hz
end
opts.baselineFrames = 1:opts.frameRate; %1s baseline. this is used for dF/F analysis later.

%% run dimensionality reduction
[bV, bU, blockInd, frameCnt, stimTime, fAlign, snap] = blockSVD(opts); %this loads raw data and does the first blockwise SVD

%% create whole-frame components
%merge dimensions if bV is in dims x trials x frames format
if iscell(bV)
    bV = cat(1,bV{:});
    if length(size(bV)) == 3
        bV = reshape(bV,size(bV,1), []);
    end
end

% combine all blue blocks and run a second SVD
[nU, s, nV] = fsvd(bV,opts.dimCnt); %combine all blocks in a second SVD
nV = s * nV'; %multiply S into V
Sv = diag(s); %keep eigenvalues

%% combine blocks back into combined components
[~, cellSize] = cellfun(@size,bU,'UniformOutput',false);
cellSize = cat(2,cellSize{:}); % get number of components in each block

% rebuild block-wise U from individual blocks
blockU = zeros(numel(snap), sum(cellSize),'single');
edgeNorm = zeros(numel(snap),1,'single');
Cnt = 0;
for iBlocks = 1 : length(bU)
    cIdx = Cnt + (1 : size(bU{iBlocks},2));
    blockU(blockInd{iBlocks}, cIdx) = blockU(blockInd{iBlocks}, cIdx) + bU{iBlocks};
    edgeNorm(blockInd{iBlocks}) = edgeNorm(blockInd{iBlocks}) + 1;
    Cnt = Cnt + size(bU{iBlocks},2);
end
edgeNorm(edgeNorm == 0) = 1; %remove zeros to avoid NaNs in blockU

%normalize blockU by dividing pixels where blocks overlap
blockU = bsxfun(@rdivide, blockU, edgeNorm);

% project block U on framewide spatial components
dSize = size(blockU);
blockU = reshape(blockU,[],dSize(end)); %make sure blockU is in pixels x componens
U = blockU * nU; %make new U with framewide components
disp('Second SVD complete'); toc;

%% do hemodynamic correction
nV = reshape(nV, size(nV,1), [], 2); % split channels
U = reshape(U,size(snap,1),size(snap,2),[]); %reshape to frame format

% do hemodynamic correction
[Vc, regC, T, hemoVar] = SvdHemoCorrect(U, nV(:,:,1), nV(:,:,2), opts.frameRate, frameCnt(2,:));

save([opts.fPath 'Vc.mat'],'Vc','U','Sv','frameCnt', 'stimTime','fAlign', '-v7.3');
save([opts.fPath 'HemoCorrection.mat'],'regC','T', 'hemoVar')
save([opts.fPath 'opts.mat'],'opts')
disp('All done!'); toc;
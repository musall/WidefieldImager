function [header,data] = loadRawData(cFile, Condition, dataType, imgSize)
% short routine to load data from WidefieldImager code.
% cFile is the file that should be opened. Condition is either
% analog to load analog data or frames to load widefield data. dataType is
% the data type. imgSize should be a variable that contains the size of the
% imaging data (as given by header = size(data)).
% Supported are tif stacks, mj2 or mp4 videos and binary files.

if ~exist('dataType','var') || isempty(dataType)
    dataType = 'uint16';
end

data = [];
header = [];
[~, ~, fileType] = fileparts(cFile); %check filetype.

if strcmpi(fileType, '.tif')
    info = imfinfo(cFile);
    data = zeros(info(1).Height, info(1).Width, info(1).SamplesPerPixel, size(info,1), 'single');
    for x = 1 : size(info,1)
        data(:,:,:,x) = imread(cFile, x);
    end
elseif strcmpi(fileType, '.mj2') || strcmpi(fileType, '.mp4')
    data = importdata(cFile);
    data = squeeze(data(:,:,1,:));
else
    fID = fopen(cFile);
    switch lower(Condition)
        
        case 'analog'
            hSize = fread(fID,1,'double'); %header size
            header = fread(fID,hSize,'double'); %Metadata. Default is: 1 = Time of Acquisition onset, 2 = Number of channels, 3 = number of values per channel
            data = fread(fID,[header(end-1),header(end)],[dataType '=>' dataType]); %get data. Last 2 header values should contain the size of the data array.
        case 'frames'
            if ~exist('imgSize','var') || isempty(imgSize) %if imgSize is not given, it should be written into the file header
                hSize = fread(fID,1,'double'); %header size
                header = fread(fID,hSize,'double'); %Metadata. Default is: 1:x = Absolute timestamps for each frame, Last 4 values: Size of each dimensions in the matrix
                imgSize = header(find(diff(header) < -1e3) + 1 : end)'; %Last 4 header values should contain the size of the data array.
            end
            data = fread(fID,[prod(imgSize),1],[dataType '=>' dataType]); %get data.
            if length(data) ~= prod(imgSize) % try again if insufficient data is found in .dat file. Sometimes fread does not get all values from file when reading from the server.
                fclose(fID);fID = fopen(cFile); %try loading data again
                hSize = fread(fID,1,'double'); %header size
                header = fread(fID,hSize,'double'); %Metadata. Defautlt is: 1:x = Absolute timestamps for each frame, Last 4 values: Size of each dimensions in the matrix
                data = fread(fID,[prod(imgSize),1],[dataType '=>' dataType]); %get data. Last 4 header values should contain the size of the data array.
            end
            data = reshape(data,imgSize); %reshape data into matrix
    end
    fclose(fID);
end
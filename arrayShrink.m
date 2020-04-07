function DataOut = arrayShrink(DataIn,mask,mode)
% Code to merge the first two dimensions of matrix 'DataIn' into one and remove
% values based on the 2D index 'mask'. The idea here is that DataIn is a stack
% of images with resolution X*Y and pixels in 'mask' should be removed to
% reduce datasize and computional load of subsequent analysis. The first 
% dimension of 'DataOut' will be the product of the X*Y of 'DataIn' minus 
% pixels in 'mask'. 
% Usage: DataOut = arrayShrink(DataIn,mask,'merge')
%
% To re-assemble the stack after computations have been done, the code
% can be called with the additional argument 'mode' set to 'split'. This
% will reconstruct the original data structure removed pixels will be
% replaced by NaNs.
% Usage: DataOut = arrayShrink(DataIn,mask,'split')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following example illustrates how to use the code when having a 3D
% image stack of size X*Y*N. X and Y are pixels and N the amount of available frames.
% Thresholding identifies the area of interest (Fig.1) and each image is compressed
% into a single vector for analzsis (Fig.2). Afterwards the original matrix
% is restored (Fig.3).
%
% Copy the example into a seperate file to test it out but do NOT uncomment
% in arrayShrink itself to ensure it's working correctly.
% 
% a = imread('coins.png'); %single image
% mask = a < mean(a(:)); %create logical index
% b = repmat(a,1,1,10); %create image stack
% b(:,:,2:2:10) = 0; %include some empty frames
% c = arrayShrink(b,mask); %merge first two dimensions
% d = arrayShrink(c,mask,'split'); %re-create original matrix b with NaNs for excluded values
% 
% figure
% subplot(1,3,1)
% imagesc(a);axis square; hold on
% contour(mask,'w','linewidth',2);
% title('Area of interest from first two dimensions')
% 
% subplot(1,3,2)
% imagesc(c'); axis square;
% set(gca,'yTick',1:2:10);set(gca,'yTickLabel',1:5);
% ylabel('Frames');xlabel('Pixels')
% title('1D vector per frame for analysis');
% 
% subplot(1,3,3)
% imagesc(d(:,:,1));axis square;
% title('Restore original matrix. Excluded data set to NaNs.')
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('mode','var')
    mode = 'merge'; %merge dimensions by default
end

dSize = size(DataIn); %size of input matrix
if dSize(1) == 1
    DataIn = squeeze(DataIn); %remove singleton dimensions
    dSize = size(DataIn);
end

if length(dSize) == 2
    if dSize(1) == 1
        DataIn = DataIn';
        dSize = size(DataIn); %size of input matrix
    end
    dSize = [dSize 1];
end

if strcmpi(mode,'merge') %merge x and y dimension into one
    
    DataIn = reshape(DataIn,[numel(mask),prod(dSize(ndims(mask)+1:end))]); %merge x and y dimension based on mask size and remaining dimensions.
    mask = mask(:); %reshape mask to vector
    DataIn(mask,:) = [];
    DataIn = reshape(DataIn,[size(DataIn,1),dSize(ndims(mask)+1:end)]);
    DataOut = DataIn;

elseif strcmpi(mode,'split') %split first dimension into x- and y- dimension based on mask size
   
    %check if datatype is single. If not will use double as a default.
    if isa(DataIn,'single')
        dType = 'single';
    else
        dType = 'double';
    end
    
    mSize = size(mask);
    mask = mask(:); %reshape mask to vector
    DataOut = NaN([numel(mask) dSize(2:end)],dType); %pre-allocate new matrix
    DataOut(~mask,:) = reshape(DataIn,sum(~mask),[]);
    DataOut = reshape(DataOut,[mSize dSize(2:end)]);
  
end
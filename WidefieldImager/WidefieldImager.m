function varargout = WidefieldImager(varargin)
% WIDEFIELDIMAGER MATLAB code for WidefieldImager.fig
%      WIDEFIELDIMAGER, by itself, creates a new WIDEFIELDIMAGER or raises the existing
%      singleton*.
%
%      H = WIDEFIELDIMAGER returns the handle to a new WIDEFIELDIMAGER or the handle to
%      the existing singleton*.
%
%      WIDEFIELDIMAGER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WIDEFIELDIMAGER.M with the given input arguments.
%
%      WIDEFIELDIMAGER('Property','Value',...) creates a new WIDEFIELDIMAGER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before WidefieldImager_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to WidefieldImager_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help WidefieldImager

% Last Modified by GUIDE v2.5 09-Aug-2020 18:16:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @WidefieldImager_OpeningFcn, ...
    'gui_OutputFcn',  @WidefieldImager_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before WidefieldImager is made visible.
function WidefieldImager_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to WidefieldImager (see VARARGIN)

% check if matlab is 2014b or newer
handles.output = hObject;
if datenum(version('-date')) < 735857 %check if matlab is 2014b or newer
    warning('Matlab version is older as 2014b. This code has not been tested on earlier versions.')
end

% some variables
handles.daqName = 'Dev3'; %name of the national instruments DAQ board
handles.extraFrames = 5; %amount of frames where the light is switched off at the end of the recording. This helps to ensure proper alignment with analog data.
handles.minSize = 100; % minimum free disk space before producing a warning (in GB)
handles.serverPath = '\\your_server_path'; %data server path

%% initialize NI card
handles = RecordMode_Callback(handles.RecordMode, [], handles); %check recording mode to create correct ni object
guidata(handles.WidefieldImager, handles);

%% initialize camera and set handles
handles.vidName = 'pcocameraadaptor'; %PCO camera
handles.vidObj = checkCamera(handles.vidName, handles); %get PCO video object

% check for other cameras if PCO adaptor is unavailable.
if isempty(handles.vidObj)
    disp('PCO camera not available. Searching for other cameras.');
    info = imaqhwinfo;
    handles.vidName = [];
    if length(info.InstalledAdaptors) == 1
        handles.vidName = info.InstalledAdaptors{1};
    elseif length(info.InstalledAdaptors) > 1 %multiple adaptors present, ask for selection
        out = listdlg('Name', 'Please select your camera', ...
            'SelectionMode','single','liststring',info.InstalledAdaptors,'listsize',[300 300]);
        if ~isempty(out)
            handles.vidName = info.InstalledAdaptors{out};
        end
    end
    handles.sBinning.Enable = 'off'; %this only works with PCO camera so disable for other cameras
end

if ~isempty(handles.vidName)
    fprintf('Using %s as current video adapter\n',handles.vidName);
    handles.vidObj = checkCamera(handles.vidName, handles); %get non-PCO video object
    preview(handles.vidObj,handles.ImagePlot.Children); %start preview
    maxRange = floor(256*0.7); %limit intensity to 70% of dynamic range to avoid ceiling effects
    cMap = gray(maxRange); cMap(end+1:256,:) = repmat([1 0 0 ],256-maxRange,1);
    colormap(handles.ImagePlot,cMap);
else
    warning('No camera found. Type "imaqhwinfo" to make sure your camera is porperly installed and running.')
end
    
if any(ismember(handles.driveSelect.String(:,1), 'g')) %start on G: drive by default
    handles.driveSelect.Value = find(ismember(handles.driveSelect.String(:,1), 'g'));
end
CheckPath(handles); %Check for data path, reset date and trialcount

% UIWAIT makes WidefieldImager wait for user response (see UIRESUME)
% uiwait(handles.WidefieldImager);


% --- Outputs from this function are returned to the command line.
function varargout = WidefieldImager_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in WaitForTrigger.
function WaitForTrigger_Callback(hObject, eventdata, handles)
% hObject    handle to WaitForTrigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of WaitForTrigger
if get(hObject, 'Value')
    set(hObject, 'String' , 'Wait for Trigger ON')
    set(hObject, 'BackgroundColor' , '[0 1 0]')
else
    set(hObject, 'String' , 'Wait for Trigger OFF')
    set(hObject, 'BackgroundColor' , '[1 0 0]')
end

% change status indicator
if handles.WaitForTrigger.Value && handles.SnapshotTaken %Waiting for trigger and snapshot taken
    handles.AcqusitionStatus.BackgroundColor = [0 1 0];
    handles.AcqusitionStatus.String = 'Waiting';
    handles.CurrentStatus.String = 'Waiting for trigger';
    AcquireData(handles.WidefieldImager); %set to acquisition mode
elseif handles.WaitForTrigger.Value && ~handles.SnapshotTaken %Waiting for trigger but no snapshot taken
    handles.AcqusitionStatus.Value = false;
    handles.AcqusitionStatus.BackgroundColor = [1 0 0];
    set(handles.AcqusitionStatus, 'String' , 'Inactive')
    handles.CurrentStatus.String = 'No snapshot taken';
elseif ~handles.WaitForTrigger.Value && handles.SnapshotTaken %Not waiting for trigger but snapshot is taken
    handles.AcqusitionStatus.Value = false;
    handles.AcqusitionStatus.BackgroundColor = [1 0 0];
    set(handles.AcqusitionStatus, 'String' , 'Inactive')
    handles.CurrentStatus.String = 'Snapshot taken';
elseif ~handles.WaitForTrigger.Value && ~handles.SnapshotTaken %Not waiting for trigger and no snapshot is taken
    handles.AcqusitionStatus.Value = false;
    handles.AcqusitionStatus.BackgroundColor = [1 0 0];
    set(handles.AcqusitionStatus, 'String' , 'Inactive')
    handles.CurrentStatus.String = 'Not ready';
end


function BaselineFrames_Callback(hObject, eventdata, handles)
% hObject    handle to BaselineFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BaselineFrames as text
%        str2double(get(hObject,'String')) returns contents of BaselineFrames as a double


% --- Executes during object creation, after setting all properties.
function BaselineFrames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BaselineFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ChangeDataPath.
function ChangeDataPath_Callback(hObject, eventdata, handles)
% hObject    handle to ChangeDataPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DataPath.String = uigetdir; %overwrites the complete file path for data storage with user selection

function DataPath_Callback(hObject, eventdata, handles)
% hObject    handle to DataPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DataPath as text
%        str2double(get(hObject,'String')) returns contents of DataPath as a double


% --- Executes during object creation, after setting all properties.
function DataPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DataPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CurrentStatus.
function CurrentStatus_Callback(hObject, eventdata, handles)
% hObject    handle to CurrentStatus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CurrentStatus


function RecordingID_Callback(hObject, eventdata, handles)
% hObject    handle to RecordingID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RecordingID as text
%        str2double(get(hObject,'String')) returns contents of RecordingID as a double


% --- Executes during object creation, after setting all properties.
function RecordingID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RecordingID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in AnimalID.
function AnimalID_Callback(hObject, eventdata, handles)
% hObject    handle to AnimalID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns AnimalID contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AnimalID

CheckPath(handles);


% --- Executes during object creation, after setting all properties.
function AnimalID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AnimalID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ExperimentType.
function ExperimentType_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ExperimentType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ExperimentType

CheckPath(handles);

% --- Executes during object creation, after setting all properties.
function ExperimentType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ExperimentType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in BlueLight.
function BlueLight_Callback(hObject, eventdata, handles)
% hObject    handle to BlueLight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.dNIdevice)
    disp('LED control not available - NI device is missing')
    set(hObject, 'Value',false)
else
    if hObject.Value
        out = false(1,3); out(handles.lightMode.Value) = true; %indicator for NI channel
        outputSingleScan(handles.dNIdevice,out)
        if handles.lightMode.Value == 1
            hObject.BackgroundColor = [0 0 1];
            hObject.String = 'BLUE is ON';
        elseif handles.lightMode.Value == 2
            hObject.BackgroundColor = [.5 0 .5];
            hObject.String = 'VIOLET is ON';
        elseif handles.lightMode.Value == 3
            hObject.BackgroundColor = [.25 0 .75];
            hObject.String = 'MIXED stim is ON';
        end
    else
        outputSingleScan(handles.dNIdevice,false(1,3))
        hObject.BackgroundColor = zeros(1,3);
        hObject.String = 'LED OFF';
    end
end

% --- Executes on key press with focus on BlueLight and none of its controls.
function BlueLight_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to BlueLight (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in TakeSnapshot.
function TakeSnapshot_Callback(hObject, eventdata, handles)
% hObject    handle to TakeSnapshot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.vidObj == 0
    disp('Snapshot not available. Check if camera is connected and restart.')
else
    h = figure('Toolbar','none','Menubar','none','NumberTitle','off','Name','Snapshot'); %create figure to show snapshot
    snap = getsnapshot(handles.vidObj); %get snapshot from video object
    temp = ls(handles.path.base); %check if earlier snapshots exist
    temp(1,8)=' '; % make sure temp has enough characters
    temp = temp(sum(ismember(temp(:,1:8),'Snapshot'),2)==8,:); %only keep snapshot filenames
    temp(~ismember(temp,'0123456789')) = ' '; %replace non-integer characters with blanks
    cNr = max(str2num(temp)); %get highest snapshot nr
    cNr(isempty(cNr)) = 0; %replace empty with 0 if no previous snapshot existed
    save([handles.path.base 'Snapshot_' num2str(cNr+1) '.mat'],'snap') %save snapshot
    imwrite(snap,[handles.path.base 'Snapshot_' num2str(cNr+1) '.jpg'], 'Bitdepth', 12) %save snapshot as jpg
    
    %     imshow(snap,'XData',[0 1],'YData',[0 1]); colormap gray; axis image;
    imshow(snap); axis image; title(['Saved as Snapshot ' num2str(cNr+1)]);
    uicontrol('String','Close','Callback','close(gcf)','units','normalized','position',[0 0 0.15 0.07]); %close button
    handles.SnapshotTaken = true; %update snapshot flag
    
    % change status indicator
    if handles.WaitForTrigger.Value %Waiting for trigger and snapshot taken
        handles.AcqusitionStatus.BackgroundColor = [0 1 0];
        handles.AcqusitionStatus.String = 'Waiting';
        handles.CurrentStatus.String = 'Waiting for trigger';
        AcquireData(handles.WidefieldImager); %set to acquisition mode
        
    elseif ~handles.WaitForTrigger.Value %Not waiting for trigger but snapshot is taken
        handles.AcqusitionStatus.Value = false;
        handles.AcqusitionStatus.BackgroundColor = [1 0 0];
        set(handles.AcqusitionStatus, 'String' , 'Inactive')
        handles.CurrentStatus.String = 'Snapshot taken';
    end
    guidata(hObject,handles);
end

% --- Executes on button press in StartPreview.
function StartPreview_Callback(hObject, eventdata, handles)
% hObject    handle to StartPreview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.vidObj == 0
    disp('Preview not available. Check if camera is connected and restart.')
else
    preview(handles.vidObj,handles.ImagePlot.Children); %start preview
    maxRange = floor(256*0.7); %limit intensity to 70% of dynamic range to avoid ceiling effects
    cMap = gray(maxRange); cMap(end+1:256,:) = repmat([1 0 0 ],256-maxRange,1);
    colormap(handles.ImagePlot,cMap);
end

% --- Executes on button press in StopPreview.
function StopPreview_Callback(hObject, eventdata, handles)
% hObject    handle to StopPreview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~(handles.vidObj == 0)
    stoppreview(handles.vidObj); %stop preview
end

% --- Executes on button press in SelectROI.
function SelectROI_Callback(hObject, eventdata, handles)
% hObject    handle to SelectROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.vidObj == 0
    disp('ROI change not available. Check if camera is connected and restart.')
else
    %% select ROI
    closepreview(handles.vidObj) %stop camera, and set new ROI position
    stop(handles.vidObj) %stop camera
    snap = getsnapshot(handles.vidObj);
    imshow(snap,'Parent',handles.ImagePlot);
    [~,ROI] = imcrop(handles.ImagePlot);
    ROI = floor(ROI);
    
    if ~isempty(ROI)
        handles.CurrentResolution.String = [num2str(ROI(3)) ' x ' num2str(ROI(4))]; %update current resolution indicator
        h = figure; %create temporary figure for saving the roi selection
%         imshow(snap,'XData',[0 1],'YData',[0 1]); hold on %plot current view
        imshow(snap); hold on %plot current view
        xVals = [repmat(ROI(1),1,2) repmat(ROI(1)+ROI(3),1,2) ROI(1)];
        yVals = [ROI(2) repmat(ROI(2)+ROI(4),1,2) repmat(ROI(2),1,2)];
        plot(xVals,yVals,'r','linewidth',2) %plot ROI outline
        savefig(gcf,[handles.path.base 'ROI.fig']) %save ROI figure
        saveas(h,[handles.path.base 'ROI.jpg']) %save ROI as jpg
        close
        
        %update preview
        handles.ROIposition = ROI;
        set(handles.vidObj,'ROIposition',handles.ROIposition);
        snap = getsnapshot(handles.vidObj); hold(handles.ImagePlot,'off');
        %         imshow(snap,[],'parent',handles.ImagePlot,'XData',[0 1],'YData',[0 1]);
        imshow(snap,[],'parent',handles.ImagePlot);
        maxRange = floor(256*0.7); %limit intensity to 70% of dynamic range to avoid ceiling effects
        cMap = gray(maxRange); cMap(end+1:256,:) = repmat([1 0 0 ],256-maxRange,1);
        colormap(handles.ImagePlot,cMap);
        preview(handles.vidObj,handles.ImagePlot.Children); %resume preview
    end
end

% --- Executes on button press in ResetROI.
function ResetROI_Callback(hObject, eventdata, handles)
% hObject    handle to ResetROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.vidObj == 0
    disp('ROI change not available. Check if camera is connected and restart.')
else
    vidRes = get(handles.vidObj,'VideoResolution');
    handles.ROIposition = [0 0 vidRes];  %default ROIposition
    stoppreview(handles.vidObj) %stop camera
    stop(handles.vidObj) %stop camera
    set(handles.vidObj,'ROIposition',handles.ROIposition); % set new ROI position
    snap = getsnapshot(handles.vidObj); %get snapshot from video object
%     imshow(snap,[],'parent',handles.ImagePlot,'XData',[0 1],'YData',[0 1]); hold(handles.ImagePlot,'off'); %plot current view
    imshow(snap,[],'parent',handles.ImagePlot); hold(handles.ImagePlot,'off'); %plot current view
    colormap(handles.ImagePlot,'gray');
    preview(handles.vidObj,handles.ImagePlot.Children); %resume preview
    maxRange = floor(256*0.7); %limit intensity to 70% of dynamic range to avoid ceiling effects
    cMap = gray(maxRange); cMap(end+1:256,:) = repmat([1 0 0 ],256-maxRange,1);
    colormap(handles.ImagePlot,cMap);
    
    % update current resolution in GUI
    imHeight = handles.ROIposition(3)-handles.ROIposition(1);
    imWidth = handles.ROIposition(4)-handles.ROIposition(2);
    handles.CurrentResolution.String = [num2str(imWidth) ' x ' num2str(imHeight)]; %update current resolution indicator
end

function CurrentResolution_Callback(hObject, eventdata, handles)
% hObject    handle to CurrentResolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CurrentResolution as text
%        str2double(get(hObject,'String')) returns contents of CurrentResolution as a double

if handles.vidObj == 0
    disp('ROI change not available. Check if camera is connected and restart.')
else
    S = get(hObject,'String'); %get resolution string
    S(strfind(S,'x'):strfind(S,'x')+1)=[]; %remove 'x' from string
    nRes = str2num(S);
    vidRes = get(handles.vidObj,'VideoResolution');
    
    if length(nRes) == 2 && nRes(1)<=vidRes(2) && nRes(2)<=vidRes(1)
        handles.ROIposition =[0 0 nRes(1) nRes(2)];
        stop(handles.vidObj) %stop camera
        set(handles.vidObj,'ROIposition',handles.ROIposition); % set new ROI position
        colormap(handles.ImagePlot,'gray');
        preview(handles.vidObj)  %resume camera preview
        maxRange = floor(256*0.7); %limit intensity to 70% of dynamic range to avoid ceiling effects
        cMap = gray(maxRange); cMap(end+1:256,:) = repmat([1 0 0 ],256-maxRange,1);
        colormap(handles.ImagePlot,cMap);
    else
        disp([get(hObject,'String') ' is not a valid input to change the resolution'])
    end
end


% --- Executes during object creation, after setting all properties.
function CurrentResolution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurrentResolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AcqusitionStatus.
function AcqusitionStatus_Callback(hObject, eventdata, handles)
% hObject    handle to AcqusitionStatus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AcqusitionStatus
if ~(get(handles.WaitForTrigger, 'value') && handles.SnapshotTaken)
    set(handles.AcqusitionStatus, 'value',false)
end

% --- Executes when user attempts to close IntrinsicImagerGUI.
function WidefieldImager_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to WidefieldImagerGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% move saved files to data save location
temp = ls(handles.path.base);

sFiles = ~cellfun('isempty',((strfind(cellstr(ls(handles.path.base)),'Snapshot_')))); %identify snapshot files
rFiles = ~cellfun('isempty',((strfind(cellstr(ls(handles.path.base)),'ROI.')))); %identify ROI file
aFiles = [temp(sFiles,:);temp(rFiles,:)]; %all files that should be moved

if ~isdir(get(handles.DataPath,'String')) && ~isempty(aFiles) %create data path if not existent and if there is data to be moved
    mkdir(get(handles.DataPath,'String'))
end
for iFiles = 1:size(aFiles,1)
    movefile([handles.path.base aFiles(iFiles,:)],[get(handles.DataPath,'String') '\' aFiles(iFiles,:)]); %move files
end

%% clear running objects
imaqreset
guidata(hObject,handles);
delete(hObject)

function AcquireData(hObject)
% function to acquire images from the camera when GUI is set to wait for
% trigger and a snapshot has been acquired
handles = guidata(hObject); %get handles

%% check for errors in trigger lines
data = false(1,3);
if ~isempty(handles.dNIdevice)
    data = inputSingleScan(handles.dNIdevice); %trigger lines
end
if sum(data) == length(data) %all triggers are active - leave acqusition mode and reset indicators
    disp('Acquire and stimulus triggers are both high - check if triggers are correctly connected and try again.')
    set(handles.WaitForTrigger, 'value' , false)
    set(handles.WaitForTrigger, 'String' , 'Wait for Trigger OFF')
    set(handles.WaitForTrigger, 'BackgroundColor' , '[1 0 0]')
    set(handles.AcqusitionStatus, 'value' , false)
    set(handles.AcqusitionStatus, 'String' , 'Inactive')
    set(handles.AcqusitionStatus, 'BackgroundColor' , '[1 0 0]')
    set(handles.AcqusitionStatus, 'String' , 'Inactive')
    handles.CurrentStatus.String = 'Snapshot taken';
    return
else
    src = getselectedsource(handles.vidObj);
    if strcmpi(handles.vidName,'pcovideoadapter')
        src.E2ExposureTime = 1000/str2double(handles.FrameRate.String) * 1000; %make sure current framerate is used
        if str2double(handles.FrameRate.String) > 10 && strcmp(handles.sBinning.String(handles.sBinning.Value),'1')
            warning('FrameRate is above 10Hz at full resolution. This can lead to performance issues.')
        end
    end
    handles.BlueLight.Value = false; BlueLight_Callback(handles.BlueLight, [], handles) %switch LED off
    colormap(handles.ImagePlot,'jet');
    
    %%move saved files to actual data save location
    temp = ls(handles.path.base);
    sFiles = ~cellfun('isempty',((strfind(cellstr(ls(handles.path.base)),'Snapshot')))); %identify snapshot files
    rFiles = ~cellfun('isempty',((strfind(cellstr(ls(handles.path.base)),'ROI.fig')))); %identify ROI file
    aFiles = [temp(sFiles,:);temp(rFiles,:)]; %all files that should be moved
    
    if ~isdir(get(handles.DataPath,'String')) %create data path if not existent
        mkdir(get(handles.DataPath,'String'))
    end
    for iFiles = 1:size(aFiles,1)
        movefile([handles.path.base aFiles(iFiles,:)],[get(handles.DataPath,'String') '\' aFiles(iFiles,:)]); %move files
    end
    save([handles.DataPath.String '\' 'handles.mat'],'handles') %save recorder handles to be able to know all the settings.
    
    % check if server location is available and create folder for behavioral data
    % 'open' indicates that this is the folder that relates to the current imaging session
    if exist(handles.serverPath,'dir')
        fPath = get(handles.DataPath,'string');
        fPath = strrep(fPath,fPath(1:2),handles.serverPath); %replace path of local drive with network drive
        fPath = [fPath '_open']; %add identifier that this session is currently being acquired (useful to point apps to the same folder)
        mkdir(fPath); %this is the server folder where other programs can write behavioral data.
    end
    
    %% start data acquisition
    stoppreview(handles.vidObj); %stop preview
    stop(handles.vidObj);
    flushdata(handles.vidObj);
    start(handles.vidObj); %get camera ready to be triggered
    StateCheck = true; %flag to control acquisition mode
    handles.lockGUI.Value = true; handles = lockGUI_Callback(handles.lockGUI, [], handles); %run callback for lock button
    
    % inactivate some parts of the GUI so they dont mess with the recording
    set(findall(handles.ExperimentID, '-property', 'enable'), 'enable', 'off')
    set(findall(handles.ControlPanel, '-property', 'enable'), 'enable', 'inactive')
    handles.FrameRate.Enable = 'off';
    handles.sBinning.Enable = 'off';
    handles.driveSelect.Enable = 'off';
    handles.ChangeDataPath.Enable = 'off';
    
    % automatically unlock software triggers if NI device is missing
    handles.triggerLock.Enable = 'on';
    if isempty(handles.dNIdevice)
        handles.triggerLock.Value = false; 
        handles = triggerLock_Callback(handles.triggerLock, [], handles); %unlock software triggers
    end
    
    while StateCheck
        %% check if still recording here
        drawnow %update GUI inputs
        StateCheck = logical(get(handles.WaitForTrigger, 'value')); %check WaitForTrigger to determine whether acqusition should still be active
        
        if ~StateCheck %leave acqusition mode and reset indicators
            disp('Acqusition stopped');
            if exist('fPath','var')
                movefile(fPath,strrep(fPath, '_open', '')); %rename data folder on the server. Later, imaging data can be moved there using 'Widefield_MoveData'.
            end
            
            % lock software trigger buttons again
            handles.triggerLock.Value = true;
            handles = triggerLock_Callback(handles.triggerLock, [], handles); %unlock software triggers
            handles.triggerLock.Enable = 'off';

            set(handles.WaitForTrigger, 'String' , 'Wait for Trigger OFF')
            set(handles.WaitForTrigger, 'BackgroundColor' , '[1 0 0]')
            set(handles.AcqusitionStatus, 'value' , false)
            set(handles.AcqusitionStatus, 'String' , 'Inactive')
            set(handles.AcqusitionStatus, 'BackgroundColor' , '[1 0 0]')
            handles.CurrentStatus.String = 'Snapshot taken';
            set(findall(handles.ExperimentID, '-property', 'enable'), 'enable', 'on')
            set(findall(handles.ControlPanel, '-property', 'enable'), 'enable', 'on')
            handles.FrameRate.Enable = 'on';
            if strcmpi(handles.vidName, 'pcocameraadaptor')
                handles.sBinning.Enable = 'on';
            end
            handles.driveSelect.Enable = 'on';
            handles.ChangeDataPath.Enable = 'on';
            CheckPath(handles);
            
            % stop camera
            stop(handles.vidObj);
            flushdata(handles.vidObj);
        
            return
        end
        
        %% wait for trial trigger
        data = false(1,3);
        if ~isempty(handles.dNIdevice)
            data = inputSingleScan(handles.dNIdevice); %trigger lines
        end
        MaxWait = str2double(get(handles.WaitingTime,'String')); %get maximum waiting time.
        aTrigger = logical(data(1)); %trigger to start data acqusition
        if sum(data) == length(data)
            aTrigger = false; %do not initate if all triggers are high together (rarely happens due to noise in the recording).
        end
        
        %check trial trigger button
        if ~aTrigger
            aTrigger = handles.TrialTrigger.Value;
            handles.TrialTrigger.Value = false;
        end
        
        if aTrigger
            set(handles.TrialNr,'String',num2str(str2double(get(handles.TrialNr,'String'))+1)); %increase TrialNr;
            
            set(handles.CurrentStatus,'String','Recording baseline'); %update status indicator
            set(handles.AcqusitionStatus, 'value', true); %set acquisition status to active
            set(handles.AcqusitionStatus, 'String' , 'Recording')
            
            if ~isempty(handles.dNIdevice)
                aID = fopen([get(handles.DataPath,'String') '\Analog_' get(handles.TrialNr,'String') '.dat'], 'wb'); %open binary file for analog data
                handles.aListen = addlistener(handles.aNIdevice,'DataAvailable', @(src, event)logAnalogData(src,event,aID,handles.AcqusitionStatus)); %listener to stream analog data to disc
                handles.aNIdevice.startBackground(); %start analog data streaming
                pause(0.2);
            end
            
            if strcmpi(handles.vidName, 'pcocameraadaptor')
                frameRate = str2double(handles.FrameRate.String);
            else
                frameRate = str2double(handles.FrameRate.String{handles.FrameRate.Value});
            end
            bSize = str2double(handles.BaselineFrames.String)*frameRate; %number of frames in baseline
            sSize = str2double(handles.PostStimFrames.String)*frameRate; %number of frames after stimulus trigger

            handles.vidObj.FramesPerTrigger = Inf; %acquire until stoppped
            trigger(handles.vidObj); %start image acquisition
            handles.BlueLight.Value = true; BlueLight_Callback(handles.BlueLight, [], handles) %switch LED on
            drawnow;
            
            tic; %timer to abort acquisition if stimulus is not received within a certain time limit
            while handles.AcqusitionStatus.Value %keep running until poststim data is recorded
                if ~isempty(handles.dNIdevice)
                    data = inputSingleScan(handles.dNIdevice); %trigger lines
                else
                    data = false(1,3);
                end
                if sum(data) == length(data)
                    stimTrigger = false; %do not proceed if all triggers are high
                else
                    stimTrigger = logical(data(2)); %trigger that indicates start of stimulus presentation
                    stopTrigger = logical(data(3)); %trigger that indicates end of trial. aborts frame acquisition of read before all poststim frames have been collected.
                    
                    %check buttons
                    if ~stimTrigger
                        drawnow;
                        stimTrigger = handles.StimTrigger.Value;
                        handles.StimTrigger.Value = false;
                    end
                    if ~stopTrigger
                        drawnow;
                        stopTrigger = handles.StopTrigger.Value;
                        handles.StopTrigger.Value = false;
                    end
                end
                
                if ~stimTrigger && (toc < MaxWait) && ~stopTrigger %record baseline frames until stimulus trigger occurs or maximum waiting time is reached
                    if handles.vidObj.FramesAvailable > bSize*2
                        getdata(handles.vidObj,bSize); %remove unnecessary frames from video capture stream
                        disp(['Waiting... Removing ' num2str(bSize) ' frames from buffer']);
                    end
                else
                    bIdx = handles.vidObj.FramesAvailable; %check available baseline frames
                    if toc < MaxWait && ~stopTrigger %record poststim if stimulus trigger was received
                        set(handles.CurrentStatus,'String','Recording PostStim');drawnow;
                        FrameWait = true;
                        while FrameWait %wait until post-stim frames are captured
                            FrameWait = handles.vidObj.FramesAvailable < (bIdx+sSize); %stop condition
                            if ~isempty(handles.dNIdevice)
                                data = inputSingleScan(handles.dNIdevice); %hardware trigger lines
                                stopTrigger = logical(data(3));
                            else
                                data = false(1,3);
                                stopTrigger = false;
                            end
                            
                            %check stop button
                            if ~stopTrigger
                                drawnow;
                                stopTrigger = handles.StopTrigger.Value;
                                handles.StopTrigger.Value = false;
                            end
                            
                            if sum(data) ~= length(data)
                                if stopTrigger %trigger that indicates if end of trial has been reached
                                    FrameWait = false;
                                    disp(['Received stop trigger. Stopped after ' num2str(handles.vidObj.FramesAvailable-bIdx) ' poststim frames']);
                                end
                            end
                        end
                    elseif stopTrigger
                        disp('Received stop trigger. No poststim data recorded'); drawnow;
                    else
                        disp('Maximum waiting time reached. No poststim data recorded'); drawnow;
                    end
                    
                    % switch off LED and grab some extra dark frames to ensure that blue and violet channels can be separated correctly.
                    recPause = (handles.extraFrames + 1) / str2double(handles.FrameRate.String); %pause long enough to get extra frames.
                    if recPause < 0.2; recPause = 0.2; end %at least 200ms pause
                    
                    handles.BlueLight.Value = false; BlueLight_Callback(handles.BlueLight, [], handles) %switch LED off
                    drawnow;
                    
                    tic
                    if ~isempty(handles.dNIdevice)
                        pause(recPause); handles.aNIdevice.stop(); %pause to ensure all analog data is written, then stop analog object
                        fclose(aID); %close analog data file
                        delete(handles.aListen); %delete listener for analog data recording
                    end
                    stop(handles.vidObj); %stop video capture
                    
                    if (bIdx-bSize) > 0
                        getdata(handles.vidObj,bIdx-bSize); %remove unnecessary frames from video object
                        disp(['Trial finished... Removing ' num2str(bIdx-bSize) ' frames from buffer']);
                    end
                    
                    if handles.vidObj.FramesAvailable < (bSize + sSize + handles.extraFrames)
                        [Data,~,frameTimes] = getdata(handles.vidObj, handles.vidObj.FramesAvailable); %collect available video data
                    else
                        [Data,~,frameTimes] = getdata(handles.vidObj,bSize + sSize + handles.extraFrames); %collect requested video data
                    end
                    frameTimes = datenum(cat(1,frameTimes(:).AbsTime)); %collect absolute timestamps
                    
                    if bIdx < bSize %if baseline has less frames as set in the GUI
                        disp(['Collected only ' num2str(bIdx) ' instead of ' num2str(bSize) ' frames in the baseline - check settings'])
                    else
                        bIdx = bSize;
                    end
                    set(handles.AcqusitionStatus, 'value', false); %stop recording
                    set(handles.AcqusitionStatus, 'String' , 'Waiting')
                    
                end
            end
            
            %% Save data to folder and clear
            set(handles.CurrentStatus,'String','Saving data');
            disp(['Trial ' get(handles.TrialNr,'String') '; Baseline Frames: ' num2str(bIdx) '; Poststim Frames: ' num2str(size(Data,4)-(bIdx)) '; Dark Frames: ' num2str(handles.extraFrames) '; Saving data ...'])
            
            % save frametimes and size of widefield data (this is useful to read binary data later)
            imgSize = size(Data);
            cFile = ([get(handles.DataPath,'String') '\frameTimes_' get(handles.TrialNr,'String') '.mat']);
            save(cFile,'frameTimes', 'imgSize'); %save frametimes
                
            if ~handles.saveTIF.Value %saw as raw binary (default because of high writing speed)
                sID = fopen([get(handles.DataPath,'String') '\Frames_' get(handles.TrialNr,'String') '.dat'], 'Wb'); %open binary stimulus file
                fwrite(sID,Data,'uint16'); %write iamging data as flat binary
                fclose(sID);
            
            else %save as tiff stack (writes slower but raw data is more accessible)
                cFile = ([get(handles.DataPath,'String') '\Frames_' get(handles.TrialNr,'String') '.tif']);
                for x = 1 : imgSize(end)
                    if length(imgSize) == 3
                        imwrite(Data(:, :, x), cFile, 'WriteMode', 'append', 'Compression', 'none');
                    elseif length(imgSize) == 4
                        imwrite(Data(:, :, :, x), cFile, 'WriteMode', 'append', 'Compression', 'none');
                    end
                end
            end
            
            % show average of current trial
            baselineAvg = squeeze(mean(Data(:,:,1,1:bIdx),4));
            stimAvg = squeeze(mean(Data(:,:,1,bIdx+1:end),4));
            stimAvg = (stimAvg-baselineAvg)./baselineAvg;
            imshow(mat2gray(stimAvg),'parent',handles.ImagePlot);
            colormap(handles.ImagePlot, parula(256));
            drawnow;
            clear Data frameTimes
            
            toc
            disp('==================================================');
            
            start(handles.vidObj); %get camera ready to be triggered again
            set(handles.CurrentStatus,'String','Waiting for trigger');
            
        end
    end
end



function PostStimFrames_Callback(hObject, eventdata, handles)
% hObject    handle to PostStimFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PostStimFrames as text
%        str2double(get(hObject,'String')) returns contents of PostStimFrames as a double


% --- Executes during object creation, after setting all properties.
function PostStimFrames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PostStimFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function TrialNr_Callback(hObject, eventdata, handles)
% hObject    handle to TrialNr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TrialNr as text
%        str2double(get(hObject,'String')) returns contents of TrialNr as a double


% --- Executes during object creation, after setting all properties.
function TrialNr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TrialNr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles = CheckPath(handles)

% Look for single-letter drives, starting at a: or c: as appropriate
ret = {};
for i = double('c') : double('z')
    if exist(['' i ':\'], 'dir') == 7
        ret{end+1} = [i ':']; %#ok<AGROW>
    end
end
handles.driveSelect.String = char(ret);

cPath = java.io.File(strtrim(handles.driveSelect.String(handles.driveSelect.Value,:)));
if (cPath.getFreeSpace / 2^30) < handles.minSize && length(handles.driveSelect.String) > 1
    answer = questdlg(['Only ' num2str((cPath.getFreeSpace / 2^30)) 'gb left on ' strtrim(handles.driveSelect.String(handles.driveSelect.Value,:)) filesep '. Change drive?'], ...
        'Drive select', 'Yes', 'No', 'No, stop asking me', 'Yes');
    if strcmp(answer, 'Yes')
        checker = true;
        cDrive = handles.driveSelect.Value; %keep current drive index
        while checker
            handles.driveSelect.Value = rem(handles.driveSelect.Value,length(handles.driveSelect.String))+1; %increase drive selection value by 1
            cPath = java.io.File(strtrim(handles.driveSelect.String(handles.driveSelect.Value,:)));
            if (cPath.getFreeSpace / 2^30) > handles.minSize
                disp(['Changed path to drive ' strtrim(handles.driveSelect.String(handles.driveSelect.Value,:)) filesep '. ' num2str((cPath.getFreeSpace / 2^30)) 'gb remaining.'])
                checker = false;
            elseif handles.driveSelect.Value == cDrive
                disp(['Could not find a drive with more then ' num2str(handles.minSize) 'gb of free space. Path unchanged.'])
                checker = false;
            end
        end
    elseif strcmp(answer, 'No, stop asking me')
        handles.minSize = 0; %don't check disk size anymore
    end
end

% set basepath and look for present animals, experiment types and past recordings
handles.path.base = [handles.driveSelect.String(handles.driveSelect.Value,:) filesep ...
    strtrim(handles.RecordMode.String{handles.RecordMode.Value}) filesep]; %set path of imaging code

if ~exist([handles.path.base 'Animals'],'dir') %check for animal path to save data
    mkdir([handles.path.base 'Animals']) %create folder if required
end

handles.AnimalID.String = cellstr(handles.AnimalID.String);
folders = dir([handles.path.base 'Animals']); %find animal folders
folders = folders([folders.isdir] & ~strncmpi('.', {folders.name}, 1));
checker = true;
for iAnimals = 1:size(folders,1) %skip first two entries because they contain folders '.' and '..'
    AllAnimals{iAnimals} = folders(iAnimals).name; %get animal folders
    if checker
        if strcmp(handles.AnimalID.String{handles.AnimalID.Value},folders(iAnimals).name) %check if current selected animal coincides with discovered folder
            handles.AnimalID.Value = iAnimals; %keep animal selection constant
            checker = false;
        end
    end
end

if isempty(iAnimals) %Check if any animals are found
    AllAnimals{1} = 'Dummy Subject'; %create dummy animal if nothing else is found
    mkdir([handles.path.base 'Animals\Dummy Subject']) %create folder for default experiment
end

handles.AnimalID.String = AllAnimals; %update AnimalID selection
if handles.AnimalID.Value > length(AllAnimals)
    handles.AnimalID.Value = 1; %reset indicator
end
if ~isempty(handles.AnimalID.Value)
    handles.path.AnimalID = AllAnimals{handles.AnimalID.Value}; %update path for current animal
end

folders = dir([handles.path.base 'Animals\' AllAnimals{handles.AnimalID.Value}]); %find Experiment folders
folders = folders([folders.isdir] & ~strncmpi('.', {folders.name}, 1));
for iExperiments = 1:size(folders,1) %skip first two entries because they contain folders '.' and '..'
    AllExperiments{iExperiments} = folders(iExperiments).name; %get experiment folders
end
if isempty(iExperiments) %Check if any experiments are found
    AllExperiments{1} = 'Default'; %create default experiment if nothing else is found
    mkdir([handles.path.base 'Animals\' AllAnimals{1} '\Default']) %create folder for default experiment
end

handles.ExperimentType.String = AllExperiments; %update experiment type selection
if size(AllExperiments,2) < handles.ExperimentType.Value; handles.ExperimentType.Value = 1; end
handles.path.ExpType = AllExperiments{handles.ExperimentType.Value}; %update path for current experiment type
cPath = [handles.path.base 'Animals\' AllAnimals{handles.AnimalID.Value} '\' AllExperiments{handles.ExperimentType.Value}]; %assign current path

if size(ls([cPath '\' date]),1) > 2 %check if folder for current date exist already and contains data
    Cnt = 1;
    while exist([cPath filesep date '_' num2str(Cnt)],'dir')
        Cnt = Cnt +1; %update counter until it is ensured that current experiment name is not used already
    end
    handles.path.RecordingID = [date '_' num2str(Cnt)]; %set folder for recording day as the date + neccesarry counter
else
    handles.path.RecordingID = date; %set folder for current recording to recording date
end
handles.RecordingID.String = handles.path.RecordingID; %update GUI
handles.DataPath.String = [cPath '\' handles.path.RecordingID]; %set complete file path for data storage
set(handles.TrialNr,'String','0'); %reset TrialNr
handles.SnapshotTaken = false; %flag for snapshot - has to be taken in order to start data acquisition
handles.CurrentStatus.String = 'Not ready'; %reset status indicator
guidata(handles.WidefieldImager, handles);


function WaitingTime_Callback(hObject, eventdata, handles)
% hObject    handle to WaitingTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WaitingTime as text
%        str2double(get(hObject,'String')) returns contents of WaitingTime as a double


% --- Executes during object creation, after setting all properties.
function WaitingTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WaitingTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function logAnalogData(src, evt, fid, flag)
% Add the time stamp and the data values to data. To write data sequentially,
% transpose the matrix.
% Modified to use an additional flag to stop ongoing data acquistion if
% false. Flag should be a handle to a control that contains a logical
% value.
%
% Example for addding function to a listener when running analog through StartBackground:
% handles.aListen = handles.aNIdevice.addlistener('DataAvailable', @(src, event)logAnalogData(src,event,aID,handles.AcqusitionStatus)); %listener to stream analog data to disc


if src.IsRunning %only execute while acquisition is still active
    if evt.TimeStamps(1) == 0
        fwrite(fid,3,'double'); %indicate number of single values in the header
        fwrite(fid,evt.TriggerTime,'double'); %write time of acquisition onset on first run
        fwrite(fid,size(evt.Data,2)+1,'double'); %write number of recorded analog channels + timestamps
        fwrite(fid,inf,'double'); %write number of values to read (set to inf since absolute recording duration is unknown at this point)
    end
    
    data = [evt.TimeStamps*1000, evt.Data*1000]' ; %convert time to ms and voltage to mV
    fwrite(fid,uint16(data),'uint16');
    %     plot(data(1,:),data(2:end,:))
    
    if ~logical(get(flag, 'value')) %check if acqusition is still active
        src.stop(); %stop recording
    end
end


% --- Executes on selection change in lightMode.
function lightMode_Callback(hObject, eventdata, handles)
% hObject    handle to lightMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lightMode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lightMode

BlueLight_Callback(handles.BlueLight, [], handles) %switch LED


% --- Executes during object creation, after setting all properties.
function lightMode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lightMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NewAnimal.
function NewAnimal_Callback(hObject, eventdata, handles)
% hObject    handle to NewAnimal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = CheckPath(handles); %Check for data path, reset date and trialcount
dPrompt = {'Enter animal ID'};
pName = 'New animal';
newMouse = inputdlg(dPrompt,pName,1,{'New mouse'});
if ~isempty(newMouse)
    mkdir([handles.path.base 'Animals' filesep newMouse{1}])
    
    dPrompt = {'Enter experiment ID'};
    pName = 'New experiment';
    
    newExp = inputdlg(dPrompt,pName,1,{[strtrim(handles.RecordMode.String{handles.RecordMode.Value}) 'Paradigm']});
    mkdir([handles.path.base 'Animals' filesep newMouse{1} filesep newExp{1}])
    
    handles = CheckPath(handles); %Check for data path, reset date and trialcount
    handles.AnimalID.Value = find(ismember(handles.AnimalID.String,newMouse{1}));
    CheckPath(handles);
end

% --- Executes on button press in NewExperiment.
function NewExperiment_Callback(hObject, eventdata, handles)
% hObject    handle to NewExperiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dPrompt = {'Enter experiment ID'};
pName = 'New experiment';
newExp = inputdlg(dPrompt,pName,1,{'New experiment'});
if ~isempty(newExp)
    mkdir([handles.path.base 'Animals' filesep handles.AnimalID.String{handles.AnimalID.Value} filesep newExp{1}])
    
    handles = CheckPath(handles); %Check for data path, reset date and trialcount
    handles.ExperimentType.Value = find(ismember(handles.ExperimentType.String,newExp{1}));
    CheckPath(handles);
end

function FrameRate_Callback(hObject, eventdata, handles)
% hObject    handle to FrameRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrameRate as text
%        str2double(get(hObject,'String')) returns contents of FrameRate as a double

if ~isempty(handles.vidObj)
    stop(handles.vidObj);
    flushdata(handles.vidObj);
    src = getselectedsource(handles.vidObj);
    if strcmpi(handles.vidName, 'pcocameraadaptor')
        if str2double(handles.FrameRate.String) > 10 && strcmp(handles.sBinning.String(handles.sBinning.Value),'1')
            answer = questdlg('Spatial binning is set to 1. This could produce a lot of data. Proceed?');
            if strcmpi(answer,'Yes')
                src.E2ExposureTime = 1000/str2double(handles.FrameRate.String) * 1000; %set current framerate
            end
        else
            src.E2ExposureTime = 1000/str2double(handles.FrameRate.String) * 1000; %set current framerate
        end
    else
        % Adjust frame rate for non-PCO camera. 
        % !! Warning !! This does not take effect with all imaq video adaptors. 
        % Make sure that your camera allows the Matlab adaptor to control the framerate.
        src.FrameRate = hObject.String{hObject.Value};
    end
end

% --- Executes during object creation, after setting all properties.
function FrameRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrameRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in sBinning.
function sBinning_Callback(hObject, eventdata, handles)
% hObject    handle to sBinning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sBinning contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sBinning

if strcmpi(handles.vidName, 'pcocameraadaptor') 
    src = getselectedsource(handles.vidObj);
    if str2double(handles.FrameRate.String) > 10 && strcmp(hObject.String(hObject.Value),'1')
        src.E2ExposureTime = 100000; %limit framerate to 10Hz
        disp('FrameRate is limited to 10Hz without spatial binning.')
    else
        src.E2ExposureTime = 1000/str2double(handles.FrameRate.String)*1000; %make sure current framerate is used
    end
    
    src.B1BinningHorizontal = hObject.String(hObject.Value);
    src.B2BinningVertical = hObject.String(hObject.Value);
    
    vidRes = get(handles.vidObj,'VideoResolution');
    handles.CurrentResolution.String = [num2str(vidRes(1)) ' x ' num2str(vidRes(2))]; %update current resolution indicator
end

% --- Executes during object creation, after setting all properties.
function sBinning_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sBinning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in RecordMode.
function handles = RecordMode_Callback(hObject, eventdata, handles)
% hObject    handle to RecordMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns RecordMode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from RecordMode

handles.ExperimentType.Value = 1;
handles.AnimalID.Value = 1;

if hObject.Value == 1 %set standard settings for widefield mapping
    handles.BaselineFrames.String = '2';
    handles.PostStimFrames.String = '23';
elseif  hObject.Value == 2 %set standard settings for behavioral recording
    handles.BaselineFrames.String = '3';
    handles.PostStimFrames.String = '6';
end
handles = CheckPath(handles);

% check NI card
daqs = daq.getDevices;
if isempty(daqs)
    disp('No NI devices found. Control acquisition with manual triggers instead.')
    handles.dNIdevice = [];
else
    if isfield(handles,'dNIdevice')
        delete(handles.dNIdevice);
        delete(handles.aNIdevice);
    end
    
    checker = false; %check if daq device with correct ID is present. Default is 'Dev1'
    for x = 1 : length(daqs)
        if strcmpi(daqs(x).ID,handles.daqName)
            checker = true;            
        end
    end
    if ~checker
        warning(['Could not find NI-DAQ: ' handles.daqName ' - Using board ' daqs(x).ID ' instead.'])
        handles.daqName = daqs(x).ID;
    end
        
    handles.dNIdevice = daq.createSession('ni'); %object for communication with NI device - digital lines
    handles.aNIdevice = daq.createSession('ni'); %object for communication with NI device - analog lines
    handles.aNIdevice.IsContinuous = true; %set to continous acquisition
    handles.aNIdevice.Rate = 1000; %set sampling rate to 1kHz
    
    addDigitalChannel(handles.dNIdevice,handles.daqName,'port1/line0:2','OutputOnly'); %output channels for blue, violet and mixed light (1.0:blue, 1.1:violet, 1.2:mixed)
    addDigitalChannel(handles.dNIdevice,handles.daqName,'port0/line0:2','InputOnly'); %input channels to trigger data acquisition and control timing of animal behavior (0.2: trial start, 0.3: stim on)
    
    ch = addAnalogInputChannel(handles.aNIdevice,handles.daqName, 0:3, 'Voltage');
    for x = 1 : length(ch)
        ch(x).TerminalConfig = 'SingleEnded'; %use single-ended recordings. Use differential when recording low SNR signals (needs to line for +/-)
    end
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function RecordMode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RecordMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in driveSelect.
function driveSelect_Callback(hObject, eventdata, handles)
% hObject    handle to driveSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns driveSelect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from driveSelect

a = strfind(handles.DataPath.String,filesep);
cPath = [strtrim(handles.driveSelect.String(handles.driveSelect.Value,:)) fileparts(handles.DataPath.String(a(1):end))];
if ~exist(cPath,'dir')
    mkdir(cPath);
end
handles = CheckPath(handles); %Check for data path, reset date and trialcount
guidata(handles.WidefieldImager, handles);

% --- Executes during object creation, after setting all properties.
function driveSelect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to driveSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in lockGUI.
function handles = lockGUI_Callback(hObject, eventdata, handles)
% hObject    handle to lockGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if hObject.Value == 0
    handles.WaitForTrigger.Enable = 'on';
    hObject.String = 'Released';
elseif hObject.Value == 1
    handles.WaitForTrigger.Enable = 'off';
    hObject.String = 'Locked';
end

% function to check camera settings
function vidObj = checkCamera(adaptorName, handles)
% check videoadapter, set settings from the GUI and start video preview. 
% Returns an empty object if 'adaptorname' is not available.
vidObj = [];
try
    imaqreset
    vidObj = videoinput(adaptorName); %get video object
    src = getselectedsource(vidObj);
    
    if strcmpi(adaptorName, 'pcocameraadaptor') %check for PCO camera and set specific settings
        src.PCPixelclock_Hz = '286000000'; %fast scanning mode
        src.E2ExposureTime = 1000/str2double(handles.FrameRate.String) * 1000; %set framerate
        src.B1BinningHorizontal = '4';
        src.B2BinningVertical = '4';
    else
        handles.FrameRate.Style = 'popupmenu'; %change menu style to indicate available frame rates
        handles.FrameRate.String = set(src, 'FrameRate');
        src.FrameRate = handles.FrameRate.String{1}; %use highest available frame rate
    end
    
    %setup and display live video feed in preview window
    vidRes = get(vidObj,'VideoResolution');
    nbands = get(vidObj,'NumberOfBands');
    handles.CurrentResolution.String = [num2str(vidRes(1)) ' x ' num2str(vidRes(2))]; %update current resolution indicator
%     imshow(zeros(vidRes(2),vidRes(1),nbands),[],'parent',handles.ImagePlot,'XData',[0 1],'YData',[0 1]); %create image object for preview
    imshow(zeros(vidRes(2),vidRes(1),nbands),[],'parent',handles.ImagePlot); %create image object for preview
    
    %set default camera configuration
    handles.ROIposition = [0 0 vidRes];  %default ROIposition
    set(vidObj,'TriggerFrameDelay',0);
    set(vidObj,'FrameGrabInterval',1);
    set(vidObj,'TriggerRepeat',0);
    set(vidObj,'ROIposition',handles.ROIposition);
    set(vidObj,'FramesPerTrigger',Inf);
    triggerconfig(vidObj,'manual');
    
catch
    vidObj = [];
end


% --- Executes on button press in TrialTrigger.
function TrialTrigger_Callback(hObject, eventdata, handles)
% hObject    handle to TrialTrigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in StimTrigger.
function StimTrigger_Callback(hObject, eventdata, handles)
% hObject    handle to StimTrigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in StopTrigger.
function StopTrigger_Callback(hObject, eventdata, handles)
% hObject    handle to StopTrigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in triggerLock.
function handles =triggerLock_Callback (hObject, eventdata, handles)
% hObject    handle to triggerLock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of triggerLock

handles.TrialTrigger.Value = false;
handles.StimTrigger.Value = false;
handles.StopTrigger.Value = false;

if hObject.Value == 0
    handles.TrialTrigger.Enable = 'on';
    handles.StimTrigger.Enable = 'on';
    handles.StopTrigger.Enable = 'on';
    hObject.String = 'Released';
elseif hObject.Value == 1
    handles.TrialTrigger.Enable = 'off';
    handles.StimTrigger.Enable = 'off';
    handles.StopTrigger.Enable = 'off';
    hObject.String = 'Locked';
end


% --- Executes on button press in saveTIF.
function saveTIF_Callback(hObject, eventdata, handles)
% hObject    handle to saveTIF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of saveTIF

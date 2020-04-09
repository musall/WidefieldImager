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
    
%     figure(50);
%     plot(data(1,:),data(2:end,:));ylim([0 200]);

    if ~logical(get(flag, 'value')); %check if acqusition is still active
        src.stop(); %stop recording
    end
end
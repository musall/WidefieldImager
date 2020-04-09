function plotAnalogData(src, evt, flag)
%short listener function to plot analog input data. flag needs to be a
%handle to a control that contains a logical value that can be set to false
%to stop plotting/acquisition.

figure(50)
plot(evt.TimeStamps,evt.Data);

if ~logical(get(flag, 'value')); %check if acqusition is still active
    src.stop(); %stop recording
end

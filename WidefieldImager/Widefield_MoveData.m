function Widefield_MoveData(path)
% short code to move imaging data from local harddrive to server.

sPath = ['Y:\smusall' path(strfind(path,':')+1:end)]; %path to server
% sPath = ['Y:\space_managed_data' path(strfind(path,':')+1:end)]; %path to server
% sPath = ['Z:\data' path(strfind(path,':')+1:end)]; %path to server
Animals = dir([path '\Animals\']);if ~isempty(Animals);Animals([1 2],:) = [];end %get animal index
for iAnimals = 1:size(Animals,1)
    Experiments = dir([path '\Animals\' Animals(iAnimals).name]);if ~isempty(Experiments);Experiments([1 2],:) = [];end %get experiment index
    for fCheck = 1:size(Experiments,1) %check for txt files
        [~,~,c] = fileparts([Experiments(fCheck).name]);
        if strcmp(c,'.txt')
            sourcePath = [path '\Animals\' Animals(iAnimals).name filesep Experiments(fCheck).name]; % path of source file
            destinationPath = [sPath '\Animals\' Animals(iAnimals).name filesep Experiments(fCheck).name]; % path of destination file
            if ~exist([sPath '\Animals\' Animals(iAnimals).name filesep], 'dir')
                mkdir([sPath '\Animals\' Animals(iAnimals).name filesep]);
            end
            copyfile(sourcePath,destinationPath); %copy text file if present
        end
    end
    dirFlags = [Experiments.isdir];Experiments = Experiments(dirFlags);% Extract only those that are directories.
    for iExp = 1:size(Experiments,1)
        Recs = dir([path '\Animals\' Animals(iAnimals).name '\' Experiments(iExp).name]);if ~isempty(Recs);Recs([1 2],:) = [];end  %get recording index
        disp([path '\Animals\' Animals(iAnimals).name '\' Experiments(iExp).name]);
        for iRecs = 1:size(Recs,1)
             base = [Animals(iAnimals).name '\' Experiments(iExp).name '\' Recs(iRecs).name]; %basic file path
             Files = dir([path '\Animals\' base]);if ~isempty(Files);Files([1 2],:) = [];end  %get recording index
             for iFiles = 1:size(Files,1)
                 if ~isdir([sPath '\Animals\' base])
                     mkdir([sPath '\Animals\' base])
                 end
                 sourcePath = [path '\Animals\' base '\' Files(iFiles).name]; % path of source file              
                 destinationPath = [sPath '\Animals\' base '\' Files(iFiles).name]; %path of destination file
                 if strcmpi('Snapshot_1.jpg',Files(iFiles).name) && exist(destinationPath, 'file') ~= 2 %dont copy if snapshot exists on server already
                     copyfile(sourcePath,destinationPath); %leave snapshot file to make sure imager programm won't record new data in the same folder
                 elseif ~strcmpi('Snapshot_1.jpg',Files(iFiles).name)
                     movefile(sourcePath,destinationPath); %move files while preserving folder structure
                 end
             end
         end
     end
end
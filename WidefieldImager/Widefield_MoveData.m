function Widefield_MoveData(localPath)
% code to move all imaging data from local harddrive to server. leaves a
% copy of the first snapshot on the local hdd.
% 'localPath' should match the base folder where the imaging data is
% located. The code expects a subfolder 'Animals' and then for each animal
% a subfolder for each paradigm (like 'BehaveParadigm') and then multiple
% folders for each recording session. The same folder structure will be
% created in the serverPath.

serverPath = '\\grid-hs\churchland_hpc_home\smusall\'; %default path to data server
serverPath = strrep(localPath,localPath(1:2),serverPath); %replace path of local drive with network drive

Animals = dir([localPath '\Animals\']);if ~isempty(Animals);Animals([1 2],:) = [];end %get animal index
for iAnimals = 1:size(Animals,1)
    Experiments = dir([localPath '\Animals\' Animals(iAnimals).name]);if ~isempty(Experiments);Experiments([1 2],:) = [];end %get experiment index
    for fCheck = 1:size(Experiments,1) %check for txt files with notes for current animal
        [~,~,c] = fileparts([Experiments(fCheck).name]);
        if strcmp(c,'.txt')
            sourcePath = [localPath '\Animals\' Animals(iAnimals).name filesep Experiments(fCheck).name]; % path of source file
            destinationPath = [serverPath '\Animals\' Animals(iAnimals).name filesep Experiments(fCheck).name]; % path of destination file
            if ~exist([serverPath '\Animals\' Animals(iAnimals).name filesep], 'dir')
                mkdir([serverPath '\Animals\' Animals(iAnimals).name filesep]);
            end
            copyfile(sourcePath,destinationPath); %copy text file if present
        end
    end
    dirFlags = [Experiments.isdir];Experiments = Experiments(dirFlags);% Extract only those that are directories.
    for iExp = 1:size(Experiments,1)
        Recs = dir([localPath '\Animals\' Animals(iAnimals).name '\' Experiments(iExp).name]);if ~isempty(Recs);Recs([1 2],:) = [];end  %get recording index
        disp([localPath '\Animals\' Animals(iAnimals).name '\' Experiments(iExp).name]);
        for iRecs = 1:size(Recs,1)
             base = [Animals(iAnimals).name '\' Experiments(iExp).name '\' Recs(iRecs).name]; %basic file path
             Files = dir([localPath '\Animals\' base]);if ~isempty(Files);Files([1 2],:) = [];end  %get recording index
             for iFiles = 1:size(Files,1)
                 if ~isdir([serverPath '\Animals\' base])
                     mkdir([serverPath '\Animals\' base])
                 end
                 sourcePath = [localPath '\Animals\' base '\' Files(iFiles).name]; % path of source file              
                 destinationPath = [serverPath '\Animals\' base '\' Files(iFiles).name]; %path of destination file
                 if strcmpi('Snapshot_1.jpg',Files(iFiles).name) && exist(destinationPath, 'file') ~= 2 %dont copy if snapshot exists on server already
                     copyfile(sourcePath,destinationPath); %leave snapshot file to make sure imager programm won't record new data in the same folder
                 elseif ~strcmpi('Snapshot_1.jpg',Files(iFiles).name)
                     movefile(sourcePath,destinationPath); %move files while preserving folder structure
                 end
             end
         end
     end
end
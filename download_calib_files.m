% this script downloads the required BaySIC calibration files
% from the official GitHub repository:
%     https://github.com/CrystalCYFu/PyBaySIC/
% It :
% 1. checks whether 2 .txt and 2 .npy 'server_*' calibration files already exist..
% 2. if found, prompts whether to re-download or stop.. if re-downloads, it
% 3. downloads the latest .zip archive of PyBaySIC
% 4. unzips it into a temporary folder
% 5. copies all "server_*" calibration files (.txt and .npy) into the current working directory
%
% run this script once before using PIP25_forward.m or PIP25_predict.m, or
% if the source repository updates the .txt and .npy files
% latest update to PyBaySIC at time of creation: Mar 2025
% 
% Written by m. osman (mo549@cam.ac.uk), Oct 2025

clear;

% first, check for existing calibration files in the current directory
txtFiles = dir(fullfile(pwd, 'server_*.txt')); % should have one dino and one bras file
npyFiles = dir(fullfile(pwd, 'server_*.npy')); % should have one dino and one bras file
if numel(txtFiles) >= 2 && numel(npyFiles) >= 2
    msg = sprintf(['Found existing calibration files in:\n%s\n', ...
                   '(%d .txt and %d .npy files)\n\n', ...
                   'Would you like to re-download the latest files from GitHub?'], ...
                   pwd, numel(txtFiles), numel(npyFiles));
    choice = questdlg(msg, 'Calibration files found!', ...
                      'Re-download', 'Stop', 'Stop');
    switch choice
        case 'Stop'
            disp('Existing calibration files detected. No download performed.');
            return;
        case 'Re-download'
            disp('Re-downloading calibration files, please wait...');
        otherwise
            disp('Indecision, indecision! No selection made... cancelling download!');
            return;
    end
end

% second ... initiate the download if you're still here! 

repoURL = 'https://github.com/CrystalCYFu/PyBaySIC/'; % make sure this is right! 
zipURL  = [repoURL 'archive/refs/heads/main.zip'];
zipFile = fullfile(pwd, 'PyBaySIC-main.zip');
destFold = fullfile(pwd, 'PyBaySIC-main');

disp('Downloading PyBaySIC from GitHub...');

try
    websave(zipFile, zipURL);
    fprintf('Downloaded ZIP to %s\n', zipFile);
catch ME
    error('Failed to download from GitHub:\n%s', ME.message);
end

fprintf('Unzipping folder...\n');
unzip(zipFile, pwd);
delete(zipFile);

if ~isfolder(destFold)
    error('Extraction failed — directory "%s" not found.', destFold);
end

% copy the folders to the current folder...
disp('Copying calibration files ...');
fileList = dir(fullfile(destFold, '**', 'server_*.*')); % find the files with the "server_*.*" tag...
if isempty(fileList)
    warning('No calibration files found in PyBaySIC repository.');
else
    for k = 1:numel(fileList)
        src = fullfile(fileList(k).folder, fileList(k).name);
        dest = fullfile(pwd, fileList(k).name);
        copyfile(src, dest);
    end
    fprintf('Copied %d calibration files to %s\n', numel(fileList), pwd);
end

% delete
disp('Cleaning up..');
rmdir(destFold, 's');

disp('DONE!');

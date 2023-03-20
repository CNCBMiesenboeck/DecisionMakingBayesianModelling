function [accuracy, reactionTimes, difficulty] = loadRawData(folder)
% helper function to open the behavioural data files
% Copyright (C) 2023 Clifford Talbot

if nargin < 1
    folder = uigetdir([pwd, filesep, 'behavioural parameters']);
end

% find accuracy and corresponding reaction time source files within parent 
% folder
accFiles = findAllFiles('trialaccuracy.csv', folder, [], [], false, true);
rtFiles = strrep(accFiles, 'trialaccuracy', 'reactiontimes');
% read files. NB task difficulty is stored in the directory structure, so 
% need to construct difficulty array based using on number of decisions
accuracy = cell2mat(cellfun(@(x) readmatrix(x), accFiles, ...
    'UniformOutput', false));
reactionTimesCell = cellfun(@(x) readmatrix(x), rtFiles, ...
    'UniformOutput', false);
reactionTimes = cell2mat(reactionTimesCell);
nDiff = cellfun(@numel, reactionTimesCell);
difficulty = zeros(size(accuracy));
nSources = numel(accFiles);
diffiucltyCounter = 0;
for ii = 1:nSources
    % get the path element corresponding to task difficulty
    diffString = getPathElement(accFiles{ii}, 1);
    % parse it
    sepIndTemp = strfind(diffString, 'vs');
    diffTemp = str2double(diffString(1:(sepIndTemp-1))) / ...
        str2double(diffString((sepIndTemp+2):end));
    % construct the portion of the difficulty array
    difficulty((1:nDiff(ii)) + diffiucltyCounter) = ...
        repmat(diffTemp, nDiff(ii), 1);
    diffiucltyCounter = diffiucltyCounter + nDiff(ii);
end
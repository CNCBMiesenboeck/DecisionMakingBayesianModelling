function [filepaths, path] = findAllFiles(pattern, path, hWaitbar, ...
    waitbarRange, returnDirs, exactFileNameFlag, startPath)
% This function recursively finds all files that match a given pattern
% within a specified directory and its subdirectories.
% Inputs:
% - pattern: a string or cell array of strings containing patterns to match
%   (wildcards are allowed).
% - path: a string specifying the directory to search in (optional).
% - hWaitbar: a handle to a waitbar object to update (optional).
% - waitbarRange: a two-element vector specifying the minimum and maximum
%   values of the waitbar (optional).
% - returnDirs: a boolean indicating whether to include directories in the
%   output (optional, defaults to false).
% - exactFileNameFlag: a boolean indicating whether to match only exact
%   filenames (optional, defaults to false).
% - startPath: a string specifying the starting directory for the
%   uigetdir dialog (optional).
% Outputs:
% - filepaths: a cell array of strings containing the paths of all files
%   that match the pattern.
% - path: a string specifying the directory that was searched.
% Copyright (C) 2023 Clifford Talbot

filepaths = cell(0, 1);
persistent filesepchar
if isempty(filesepchar)
    filesepchar = filesep;
end
if nargin < 2 || ~ischar(path) || ~exist(path, 'dir')
    if nargin < 7
        startPath = '';
    end
    path = uigetdir(startPath);
    if length(path) == 1
        return;
    end
elseif path(end) == filesepchar
    path(end) = [];
end

if nargin >= 3 
    if ~isempty(hWaitbar) && isgraphics(hWaitbar) && ~isnumeric(hWaitbar) ...
            && isvalid(hWaitbar)
        waitbarOn = true;
        if nargin < 4
            waitbarRange = [0, 1];
        end
        waitbarMin = waitbarRange(1);
        waitbarRange = waitbarRange(2) - waitbarRange(1);
        waitbarResolution = 50;
        if waitbarRange < (1/waitbarResolution)
            waitbarOn = false;
        end
    else
        waitbarOn = false;
        waitbarMin = 0;
        waitbarRange = 1;
    end
else
    waitbarOn = false;
    waitbarMin = 0;
    waitbarRange = 1;    
    hWaitbar = gobjects; delete(hWaitbar);
end

if nargin < 5
    returnDirs = false;
end

if nargin < 6
    exactFileNameFlag = false;
end

if exactFileNameFlag
    hFunMatch = @(x, y) strcmp(x, y);
else
    hFunMatch = @contains;
end

if ~iscell(pattern)
    pattern = {pattern};
end

contents = dir(path);
contents(strcmp('.', {contents.name})...
    | strcmp('..', {contents.name})) = [];

paths = {contents.name};
isdirResult = [contents.isdir];
nPaths = length(paths);

if waitbarOn && isvalid(hWaitbar)
    waitbarUpdateRate = max(1, floor(nPaths/(waitbarRange*waitbarResolution)));
    waitbar(waitbarMin, hWaitbar, 'finding files');
end

for p = 1:nPaths
    testpath = [path, filesepchar, paths{p}];
%     isdirResult = exist(testpath, 'dir');
    if isdirResult(p) || returnDirs
        if returnDirs
            test = true;
            for pat = 1:length(pattern)
                test = test && hFunMatch(paths{p}, pattern{pat}); %~isempty(strfind(paths{p}, pattern{pat}));
            end
            if test
                filepaths{end+1, 1} = testpath;
            end
        end
        
        if isdirResult(p)
            filepaths = [filepaths; findAllFiles(pattern, testpath, ...
                hWaitbar, waitbarMin + waitbarRange*[(p-1)/nPaths, p/nPaths], ...
                returnDirs, exactFileNameFlag)];
        end
    elseif length(pattern) == 1 && isempty(pattern{1})
        filepaths{end+1, 1} = testpath;
    else
        test = true;
        for pat = 1:length(pattern)
            test = test && hFunMatch(paths{p}, pattern{pat}); % ~isempty(strfind(paths{p}, pattern{pat}));
        end
        if test
            filepaths{end+1, 1} = testpath;
        end
    end
    
    if waitbarOn 
        if isvalid(hWaitbar)
            chktemp = (p/waitbarUpdateRate);
            if floor(chktemp) == chktemp
                waitbar(waitbarMin + waitbarRange*p/nPaths, hWaitbar);
            end
        else
            break;
        end
    end
end
end
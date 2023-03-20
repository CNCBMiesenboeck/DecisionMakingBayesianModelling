function [str, upto, after] = getPathElement(path, element)
% returns a folder in the dataPath element steps up from path
% element = 0 returns the data file name
% element = 1 returns the parent folder
% etc
% str = the single element
% upto = the path upto and including the element
% Copyright (C) 2023 Clifford Talbot

ind = find(path==filesep);
if element >= numel(ind)
    str = '';
    upto = '';
    after = path;
elseif element == 0
    upto = path(1:end);
    str = upto(ind(end)+1:end);
    after = '';
else
    element = element - 1;
    upto = path(1:ind(end-element)-1);
    str = upto(ind(end-element-1)+1:end);
    after = path(ind(end-element)+1:end);
end
end
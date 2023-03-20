function [hFcnOut, initOut, hIndFcnsOut, paramsInds] = globalHandle(hFns, init, config)
% This is (sort of) a helper function that generates a handle to the target
% pdf (hFcnOut) and the starting point (initOut) of a global model.
% hFns is a cell array of the source functions i.e.  the handles to the 
% pdfs of the optostim-independent models for the data at the individual 
% stimulation levels of the genotype.
% init is an array defining the starting point for the four fundamental 
% parameters i.e. init(1) = A, (2) = k, (3) = tND, (4) = sigmaND
% config selects the global model:
%   1 = free A, k. global tND, sigmaND;     2*SS+2 independent parameters
%   2 = free k. global A, tND, sigmaND;     SS+3   independent parameters
%   3 = free A. global k, tND, sigmaND;     SS+3   independent parameters
%   4 = global A, k, tND, sigmaND;          4      independent parameters
%                                               SS = number of optostim
%                                                    intensies
% The output function (i.e. the global model) remaps the independent input 
% parameters to the source functions and sums the output (log-probability)

% empty elements of hFns arise because I started writing the code assuming 
% that there would be the same number of optostim intensities for all of 
% the genotypes. Delete the empty elements.
% Copyright (C) 2023 Clifford Talbot

invalidFns = cellfun(@isempty, hFns);
noHighIntFlag = any(invalidFns);
hFns(invalidFns) = [];

% this isn't the most versatile way to code it. But for the data set we
% have where there are either two or three stimulation levels, it is fine
% and perhaps easier to understand. It will require modification for
% alternate combinations.
switch config
    case 1 % fixed tND, sigmaND
        if ~noHighIntFlag
            p2fIndex = [1, 2, 1, 2, 1, 2, 3, 4];
            paramsInds = [...
                1, 2, 7, 8; ...
                3, 4, 7, 8; ...
                5, 6, 7, 8];
        else
            p2fIndex = [1, 2, 1, 2, 3, 4];
            paramsInds = [...
                1, 2, 5, 6; ...
                3, 4, 5, 6];
        end
    case 2 % fixed A, tND, sigmaND
        if ~noHighIntFlag
            p2fIndex = [1, 2, 2, 2, 3, 4];
            paramsInds = [...
                1, 2, 5, 6; ...
                1, 3, 5, 6; ...
                1, 4, 5, 6];
        else
            p2fIndex = [1, 2, 2, 3, 4];
            paramsInds = [...
                1, 2, 4, 5; ...
                1, 3, 4, 5];
        end
    case 3 % fixed k, tND, sigmaND
        if ~noHighIntFlag
            p2fIndex = [1, 1, 1, 2, 3, 4];
            paramsInds = [...
                1, 4, 5, 6; ...
                2, 4, 5, 6; ...
                3, 4, 5, 6];
        else
            p2fIndex = [1, 1, 2, 3, 4];
            paramsInds = [...
                1, 3, 4, 5; ...
                2, 3, 4, 5];
        end
    case 4 % fixed A, k, tND, sigmaND
        if ~noHighIntFlag
            p2fIndex = [1, 2, 3, 4];
            paramsInds = [...
                1, 2, 3, 4;
                1, 2, 3, 4;
                1, 2, 3, 4];
        else
            p2fIndex = [1, 2, 3, 4];
            paramsInds = [...
                1, 2, 3, 4; ...
                1, 2, 3, 4];
        end
end
if numel(init) > 4
    p2fIndex(end+1) = 5;
    paramsInds(:, end+1) = max(paramsInds, [], 'all')+1;
end

initOut = init(p2fIndex);
hFcnOut = @(params) sum(arrayfun(@(ii) hFns{ii}(params(paramsInds(ii, :))), ...
    1:numel(hFns)));
hIndFcnsOut = arrayfun(@(ii) @(params) hFns{ii}(params(paramsInds(ii, :))), ...
    1:numel(hFns), 'UniformOutput', false);
if noHighIntFlag
    hIndFcnsOut = [{[]}, hIndFcnsOut];
end
end
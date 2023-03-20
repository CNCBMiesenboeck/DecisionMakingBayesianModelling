% % % Copyright (C) 2023 Clifford Talbot
% % % This program is free software: you can redistribute it and/or modify
% % % it under the terms of the GNU General Public License as published by
% % % the Free Software Foundation, either version 3 of the License, or
% % % (at your option) any later version.
% % % 
% % % This program is distributed in the hope that it will be useful,
% % % but WITHOUT ANY WARRANTY; without even the implied warranty of
% % % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% % % GNU General Public License for more details.
%% Section 1 
% This section loads data and initialises arrays for processing.
% After running this section, either:
%   run a new analysis by proceeding to section 2
%   load a previously saved analysis by loading the "fitting" file and then 
%   run section 3 to generate functions to handles.

clear

% select data parent directory
dataParentDIR = uigetdir('behavioural parameters');

% data from each genotype are sepatated into the following folders
sources = {...
    'NP6024 AMPAR (targeted Chrimson)', ...
    'NP6024 CS (no Chrim)', ...
    'NP6024 NLG1 (untargeted Chrimson)', ...
    'VT030604 AMPAR (prime KC)'};
nGT = numel(sources);

% there are four global fit (holding none, either, or both of A/k (bound
% height/drift rate) constant accross the optical stimulation intensity
% levels. 
nGlobalFitOpts = 4;

% data from each opto-stim level are sepatated into the following folders
intFolders = {'high intensity', 'low intensity', 'no light'};
nInt = numel(intFolders);

% variables to save following computation
saveVars = {'sources', 'intFolders', 'rndDDMindep', 'rndEXDindep', 'dt', ...
    'initEXD', 'initDDM', 'rndEXDglobal', 'rndDDMglobal', 'globalInitEXD', ...
    'globalInitDDM', 'timeElapsedIndep', 'timeElapsedGlob', 'burninDDM', ...
    'burninEXD', 'nSamplesPP', 'nDataIndep', 'nDataGlob', 'GTRange', ...
    'intRange', 'globalFitInds', 'loadedFlag', 'totalTimeElapsed', ...
    'fitIndepFlag', 'BICEXDIndep', 'BICDDMIndep', 'BICEXDGlob', ...
    'BICDDMGlob', 'filePath', 'seedRndGenFlag', 'rndSeed'};

% create empty arrays for source data/results. The dimension with: nGT
% corresponds to the genotype, nInt = optostim intensity.
% acc = trial accuracy: 1 for correct, 0 for incorrect.
% rt = reaction time
% diff = odour contration ratio (closer to 1 = more difficult)
% logDiff = |log(diff)| (used as numerical quantification of difficulty)
[acc, rt, diff, logDiff] = deal(cell(nGT, nInt));

% rndXXXindep = output from slicesample for analysis with all parameters
% independent of optostim intensity level.
% hLogXXXindep = handles to target pdf
% (EXD = extrema detection, DDM = drift diffusion model)
[rndDDMindep, rndEXDindep, hLogEXDindep, hLogDDMindep] = ...
    deal(repmat({cell(nInt, 1)}, nGT, 1));

% rndXXXglobal, hLogXXXGlob = as above, but for global parameters
% globalInitXXX = starting point for global analyses (used for convenience)
[rndDDMglobal, rndEXDglobal, hLogEXDGlob, hLogDDMGlob, ...
    globalInitEXD, globalInitDDM] = deal(cell(nGT, nGlobalFitOpts));
[globalInitEXDSource, globalInitDDMSource] = deal(cell(nGT, 1));

% nDataIndep = number of data points (decisions) at each intensity and
% for each genotype
% nDataGlob = number of data points for the global models
% BICXXXIndep/Glob = BIC values
[nDataIndep, BICEXDIndep, BICDDMIndep] = deal(repmat({zeros(nInt, 1)}, nGT, 1));
[BICEXDGlob, BICDDMGlob] = deal(zeros(nGT, nGlobalFitOpts));
nDataGlob = zeros(nGT, 1);

% used to measure the time of computation
timeElapsedIndep = repmat({zeros(2, nInt)}, nGT, 1);
timeElapsedGlob = repmat({zeros(2, nGlobalFitOpts)}, nGT, 1);

% load the raw data
loadedFlag = false(nGT, nInt);
for gtInd = 1:nGT
    for intInd = 1:nInt
        loadSource = [dataParentDIR, filesep, sources{gtInd}, filesep, intFolders{intInd}];
        if isfolder(loadSource)
            [acc{gtInd, intInd}, rt{gtInd, intInd}, diff{gtInd, intInd}] = ...
                loadRawData(loadSource);
            logDiff{gtInd, intInd} = abs(log10(diff{gtInd, intInd}));
            loadedFlag(gtInd, intInd) = true;
            nDataIndep{gtInd}(intInd) = numel(acc{gtInd, intInd});
        end
    end
end

% flag used to signal create handles to global functions only (i.e. not
% running new analysis)
createHandlesOnly = true;

%% section 2
% run this section (modify parameters as necessary) before running section 
% 3 to generate a new sample set.

% change flag to indicate new analysis
createHandlesOnly = false;

% to fix the random number generator seed.
seedRndGenFlag = true;
rndSeed = 1; % set to 1 for the analysis used in the publication

% if the output from only the global models is required set this to false.
fitIndepFlag = true;

% set this to true to test the script
testingFlag = false;
if testingFlag
    burninDDM = 5;
    burninEXD = 5;
    nSamplesPP = 20;
else
    burninDDM = 300; % no. samples for burn-in. Can be set independently 
    burninEXD = 300; % for the two models
    nSamplesPP = 500; % no. samples per parameter for output pdf estimate
end

% the genotypes to use in the analysis (indicies of the "sources" cell
% array defined in section 1)
GTRange = 1:nGT;

% the optostim levels to use (defined in intFolders)
intRange = 1:nInt;

% to select the global models to use:
%   1 = free A, k. global tND, sigmaND
%   2 = free k. global A, tND, sigmaND
%   3 = free A. global k, tND, sigmaND
%   4 = global A, k, tND, sigmaND
globalFitInds = 1:4;

% snapshot time window for extrema detection model
dt = 0.0005;

% starting parameters for the slice sampling.
% these starting points were obtained during initial trial runs on an 
% earlier data set (a subset of the final data set). They can be changed 
% quite a bit and still give the same result (within errors), provided that
% the burn in is sufficient. But these particular values are required to 
% reproduce the sampled pdfs presented in the publication.
% element 1, 2, 3, 4: A, k, tND, sigmaND.
initEXD = [0.0777, 22, 0.25, 0.07]; % extrema detection
initDDM = [1.2, 0.94, 0.2, 0.023];  % drift diffusion

% number of samples in the sampled pdf for the independent models
nSamplesIndep = nSamplesPP*numel(initEXD); 
%% Section 3
timeTotalStart = tic;
parfor gtInd = GTRange % loop over genotypes
% for gtInd = GTRange % toggle this line with the one above to parallel computing
    % initialise random number generator if necessary
    if ~createHandlesOnly && seedRndGenFlag
        rng(rndSeed);
    end
    
    % loop over different optostim intensity levels for the independent
    % models
    for intInd = intRange
        if loadedFlag(gtInd, intInd)
%%
% I sometimes "reset" the indentation for large chunks of code at the
% innermost loop
indexString = sprintf('GT: %d, int: %d, nData: %d', gtInd, intInd, nDataIndep{gtInd}(intInd));

%%%% estimate pdf for the extrema detection model
% create handle to the model for the data subset (genotype and optostim 
% intensity) as required for the slicesample function.
% params = 4 element array corresponding to A, k, tND, sigmaND
hLogEXDindep{gtInd}{intInd} = @(params) sum(logP_RT_acc(rt{gtInd, intInd}, ...
    acc{gtInd, intInd}==1, ...
    logDiff{gtInd, intInd}, ...
    params(1), params(2), params(3), params(4), dt));
% ensure that starting point doesn't yield p = 0 (or log(p)=-inf)
testP = ~isinf(hLogEXDindep{gtInd}{intInd}(initEXD));
if ~testP
    % if invalid starting point print error message
    fprintf('indep EXD %s choose different init\n', indexString);
elseif ~createHandlesOnly && fitIndepFlag
    % otherwise run the slice sample algorithm
    fprintf('indep EXD %s running... ', indexString); % display update
    tStartTemp = tic; % start timer
    rndEXDindep{gtInd}{intInd} = slicesample(initEXD, nSamplesIndep, ...
        'logpdf', hLogEXDindep{gtInd}{intInd}, 'burnin', burninEXD);
    timeElapsedIndep{gtInd}(1, intInd) = toc(tStartTemp); % store elapsed time
    paramsEXDTemp = densityMax(rndEXDindep{gtInd}{intInd}); % compute best estimates 
    BICEXDIndep{gtInd}(intInd) = numel(paramsEXDTemp)*log(nDataIndep{gtInd}(intInd)) - ...
        2*hLogEXDindep{gtInd}{intInd}(paramsEXDTemp); % compute BIC
    % display output
    fprintf('complete in %.2fs. BIC = %.2e. estimates: ', ...
        timeElapsedIndep{gtInd}(1, intInd), BICEXDIndep{gtInd}(intInd));
    fprintf('%.2e, ', paramsEXDTemp);
    fprintf('\n');
end

%%%% independent drift diffusion model
% same procedure as above, but for DDM
hLogDDMindep{gtInd}{intInd} = @(params) sum(pDDM(rt{gtInd, intInd}, ...
    acc{gtInd, intInd}==1, ...
    logDiff{gtInd, intInd}, ...
    params(1), params(2), params(3), params(4)));
testP = ~isinf(hLogDDMindep{gtInd}{intInd}(initDDM));
if ~testP
    fprintf('indep DDM %s choose different init\n', indexString);
elseif ~createHandlesOnly && fitIndepFlag
    fprintf('indep DDM %s running... ', indexString);
    tStartTemp = tic;
    rndDDMindep{gtInd}{intInd} = slicesample(initDDM, nSamplesIndep, ...
        'logpdf', hLogDDMindep{gtInd}{intInd}, 'burnin', burninDDM);
    timeElapsedIndep{gtInd}(2, intInd) = toc(tStartTemp);
    paramsDDMTemp = densityMax(rndDDMindep{gtInd}{intInd});
    BICDDMIndep{gtInd}(intInd) = numel(paramsDDMTemp)*log(nDataIndep{gtInd}(intInd)) - ...
        2*hLogDDMindep{gtInd}{intInd}(paramsDDMTemp);
    fprintf('complete in %.2fs. BIC = %.2e. estimates: ', ...
        timeElapsedIndep{gtInd}(2, intInd), BICDDMIndep{gtInd}(intInd));
    fprintf('%.2e, ', paramsDDMTemp);
    fprintf('\n');
end
        end
    end
    % end of independent model section
    
    %%%% global models section
    % generate the starting point of the global models based on the binned 
    % posterior pdfs of the independent models
    if ~createHandlesOnly && fitIndepFlag
        % get posterior pdfs for the existing optostim levels and combine
        % them into temp. 
        temp = [rndEXDindep{gtInd}{loadedFlag(gtInd, :)}];
        % columns [1, 5, (9)], [2, 6, (10)] etc of temp contain the 
        % estimates of A, k etc at the different optostim intensities.
        % Permute and reshape temp, then find the maximum
        globalInitEXDSource{gtInd} = densityMax(reshape(...
            temp(:, [1:4:end, 2:4:end, 3:4:end, 4:4:end]), [], 4));
        % do the same for the DDM
        temp = [rndDDMindep{gtInd}{loadedFlag(gtInd, :)}];
        globalInitDDMSource{gtInd} = densityMax(reshape(...
            temp(:, [1:4:end, 2:4:end, 3:4:end, 4:4:end]), [], 4));
    else
        % the independent models were not run, then just use the starting
        % user defined starting point from section 2
        globalInitEXDSource{gtInd} = initEXD;
        globalInitDDMSource{gtInd} = initDDM;
    end
    
    % number of data points for the global models = all data points for the
    % different optostim intensities
    nDataGlob(gtInd) = sum(nDataIndep{gtInd}, 'all');
    
    % loop over the global models
    for cc = globalFitInds
        indexString = sprintf('GT: %d, glob: %d, nData: %d', gtInd, cc, nDataGlob(gtInd));
        
        % generate the handle to the target pdf of the global model and the
        % starting point (which is based on globalInitEXDSource, with
        % duplicate elements for the non-global parameters).
        [hLogEXDGlob{gtInd, cc}, globalInitEXD{gtInd, cc}] = ...
            globalHandle(hLogEXDindep{gtInd}, globalInitEXDSource{gtInd}, cc);        
        % total number of samples in sampled pdf
        nSamplesLocal = nSamplesPP * numel(globalInitEXD{gtInd, cc});
        % ensure that the starting point gives non zero probability
        testP = ~isinf(hLogEXDGlob{gtInd, cc}(globalInitEXD{gtInd, cc}));
        % if the starting point is invalid, try generating one based on the
        % user defined starting point from section 2
        if ~testP && fitIndepFlag
            [~, globalInitEXD{gtInd, cc}] = globalHandle(...
                hLogEXDindep{gtInd}, initEXD, cc);
            testP = ~isinf(hLogEXDGlob{gtInd, cc}(globalInitEXD{gtInd, cc}));
        end
        if ~testP
            % if starting point is still invalid, print error message
            fprintf('global EXD %s choose different init\n', indexString);
        elseif ~createHandlesOnly
            % otherwise run slice sample algorithm 
            fprintf('global EXD %s running... ', indexString); % display update
            tStartTemp = tic; % start timer
            rndEXDglobal{gtInd, cc} = slicesample(globalInitEXD{gtInd, cc}, ...
                nSamplesLocal, 'logpdf', hLogEXDGlob{gtInd, cc}, ...
                'burnin', burninEXD); % call slicesample
            timeElapsedGlob{gtInd}(1, cc) = toc(tStartTemp); % store elapsed time
            paramsEXDTemp = densityMax(rndEXDglobal{gtInd, cc}); % compute best estimates
            BICEXDGlob(gtInd, cc) = numel(paramsEXDTemp)*log(nDataGlob(gtInd)) - ...
                2*hLogEXDGlob{gtInd, cc}(paramsEXDTemp); % compute BIC
            % display output
            fprintf('complete in %.2fs. BIC = %.2e. estimates: ', ...
                timeElapsedGlob{gtInd}(1, cc), BICEXDGlob(gtInd, cc));
            fprintf('%.2e, ', paramsEXDTemp);
            fprintf('\n');
        end
        
        % same for DDM
        [hLogDDMGlob{gtInd, cc}, globalInitDDM{gtInd, cc}] = ...
            globalHandle(hLogDDMindep{gtInd}, globalInitDDMSource{gtInd}, cc);
        testP = ~isinf(hLogDDMGlob{gtInd, cc}(globalInitDDM{gtInd, cc}));
        if ~testP && fitIndepFlag
            [~, globalInitDDM{gtInd, cc}] = globalHandle(...
                hLogDDMindep{gtInd}, initDDM, cc);
            testP = ~isinf(hLogDDMGlob{gtInd, cc}(globalInitDDM{gtInd, cc}));
        end
        if ~testP
            fprintf('global DDM %s choose different init\n', indexString);
        elseif ~createHandlesOnly
            fprintf('global DDM %s running... ', indexString);
            tStartTemp = tic;
            rndDDMglobal{gtInd, cc} = slicesample(globalInitDDM{gtInd, cc}, ...
                nSamplesLocal, 'logpdf', hLogDDMGlob{gtInd, cc}, ...
                'burnin', burninDDM);
            timeElapsedGlob{gtInd}(2, cc) = toc(tStartTemp);
            paramsDDMTemp = densityMax(rndDDMglobal{gtInd, cc});
            BICDDMGlob(gtInd, cc) = numel(paramsDDMTemp)*log(nDataGlob(gtInd)) - ...
                2*hLogDDMGlob{gtInd, cc}(paramsDDMTemp);
            fprintf('complete in %.2fs. BIC = %.2e. estimates: ', ...
                timeElapsedGlob{gtInd}(2, cc), BICDDMGlob(gtInd, cc));
            fprintf('%.2e, ', paramsDDMTemp);
            fprintf('\n');
        end
    end
    fprintf('\n');
end

if createHandlesOnly
    % if loading a saved posterior pdf, compute the best estimates and BICs
    if (all(cellfun(@(x) all(x==0, 'all'), BICDDMIndep), 'all') && ...
            all(BICDDMGlob==0, 'all')) % || recomputeBICFlag % recomputeBIC if changing densityMax function
        for gtInd = 1:4
            for intInd = 1:3
                if loadedFlag(gtInd, intInd)
                    paramsDDM = densityMax(rndDDMindep{gtInd}{intInd}(:, :));
                    paramsEXD = densityMax(rndEXDindep{gtInd}{intInd}(:, :));
                    k = numel(initEXD)*log(nDataIndep{gtInd}(intInd));
                    BICEXDIndep{gtInd}(intInd) = k - 2*hLogEXDindep{gtInd}{intInd}(paramsEXD);
                    BICDDMIndep{gtInd}(intInd) = k - 2*hLogDDMindep{gtInd}{intInd}(paramsDDM);
                end
                for cc = 1:4
                    paramsDDM = densityMax(rndDDMglobal{gtInd, cc}(:, :));
                    paramsEXD = densityMax(rndEXDglobal{gtInd, cc}(:, :));
                    k = numel(globalInitEXD{gtInd, cc})*log(nDataGlob(gtInd));
                    BICEXDGlob(gtInd, cc) = k - 2*hLogEXDGlob{gtInd, cc}(paramsEXD);
                    BICDDMGlob(gtInd, cc) = k - 2*hLogDDMGlob{gtInd, cc}(paramsDDM);
                end
            end
        end
    end
    %% print the output from the saved file, or run this section to
    % re-output the messages
    for gtInd = GTRange
        for intInd = intRange
            if loadedFlag(gtInd, intInd)
                indexString = sprintf('GT: %d, int: %d, nData: %d', gtInd, intInd, nDataIndep{gtInd}(intInd));
                fprintf('indep EXD %s ', indexString);
            fprintf('complete in %.2fs. BIC = %.1f. estimates: ', ...
                timeElapsedIndep{gtInd}(1, intInd), BICEXDIndep{gtInd}(intInd));
                fprintf('%.2e, ', densityMax(rndEXDindep{gtInd}{intInd}));
                fprintf('\n');

                fprintf('indep DDM %s ', indexString);
                fprintf('complete in %.2fs. BIC = %.1f. estimates: ', ...
                    timeElapsedIndep{gtInd}(2, intInd), BICDDMIndep{gtInd}(intInd));
                fprintf('%.2e, ', densityMax(rndDDMindep{gtInd}{intInd}));
                fprintf('\n');
            end
        end
        for cc = globalFitInds
            indexString = sprintf('GT: %d, glob: %d, nData: %d', gtInd, cc, nDataGlob(gtInd));
            fprintf('global EXD %s ', indexString);
            fprintf('complete in %.2fs. BIC = %.1f. estimates: ', ...
                timeElapsedGlob{gtInd}(1, cc), BICEXDGlob(gtInd, cc));
            fprintf('%.2e, ', densityMax(rndEXDglobal{gtInd, cc}));
            fprintf('\n');
            fprintf('global DDM %s ', indexString);
            fprintf('complete in %.2fs. BIC = %.1f. estimates: ', ...
                timeElapsedGlob{gtInd}(2, cc), BICDDMGlob(gtInd, cc));
            fprintf('%.2e, ', densityMax(rndDDMglobal{gtInd, cc}));
            fprintf('\n');
        end
        fprintf('\n');
    end
else
    totalTimeElapsed = toc(timeTotalStart);
    % save the results and relevant parameters
    filePath = ['fitting ', datestr(datetime, 'yy-mm-dd HHMM'), '.mat'];
    save(filePath, saveVars{:});
end
fprintf('total time: %.2f, parallel time: %.2f\n', totalTimeElapsed, ...
    sum([timeElapsedIndep{:}], 'all') + sum([timeElapsedGlob{:}], 'all'));

%% Section 4. To export data to text files. 
% this saves the posterior pdf samples for the DDM model only to tab 
% delimited text files and also summary sheets (tab delimited text files)
% for both DDM and EXD models. The script creates a folder in the working 
% directory named "output {date stamp}" where date stamp is in the format 
% [YY-MM-DD hhmm]. The output files are saved in this new folder. The file 
% names should be easy to interpret.

saveDIR = [pwd, filesep, datestr(datetime, 'output yy-mm-dd HHMM'), filesep];
mkdir(saveDIR);
globalTitles = {...
    'M0. Individual: bound, drift, tND, sigmaND.', ...
    'M1. Individual: bound, drift. Global: tND, sigmaND.', ...
    'M2. Individual: drift. Global: A, tND, sigmaND.', ...
    'M3. Individual: bound. Global: drift, tND, sigmaND.', ...
    'M4. Global: bound, drift, tND, sigmaND.'};

paramsTitles = {'bound', 'drift', 'tND', 'sigmaND', 'bias', 'BIC'};
paramsIndsMaps = {...
    [1, 2, 3, 4], ...
    [1, 2, 1, 2, 1, 2, 3, 4], ...
    [1, 2, 2, 2, 3, 4], ...
    [1, 1, 1, 2, 3, 4], ...
    [1, 2, 3, 4]};

paramsIndsNoBIC = paramsIndsMaps;
paramsIndsMaps = cellfun(@(x) [x, 6], paramsIndsMaps, ...
    'UniformOutput', false);
paramsIndsMaps{1} = repmat(paramsIndsMaps{1}, 1, 3);

[paramsEXDindep, paramsDDMindep, stdEXDindep, stdDDMindep, ...
    c95EXDindep, c95DDMindep] = deal(cell(nGT, nInt));
[paramsEXDGlob, paramsDDMGlob, stdEXDGlob, stdDDMGlob, ...
    c95EXDGlob, c95DDMGlob] = deal(cell(nGT, nGlobalFitOpts));

for ss = GTRange
    for ii = intRange
        if loadedFlag(ss, ii)
            [paramsEXDindep{ss, ii}, stdEXDindep{ss, ii}, ...
                c95EXDindep{ss, ii}] = densityMax(rndEXDindep{ss}{ii});
            [paramsDDMindep{ss, ii}, stdDDMindep{ss, ii}, ...
                c95DDMindep{ss, ii}] = densityMax(rndDDMindep{ss}{ii});
            
            fid = fopen([saveDIR, 'DDM ', strrep(globalTitles{1}, ':', ''), ...
                ' ', sources{ss}, ' ', intFolders{ii}, '.txt'], 'w');
            fprintf(fid, '%s\t', paramsTitles{paramsIndsNoBIC{1}});
            fprintf(fid, '\n');
            fprintf(fid, [repmat('%.4e\t', 1, size(rndDDMindep{ss}{ii}, 2)), '\n'], ...
                rndDDMindep{ss}{ii}');
            fclose(fid);
        end
    end
    for cc = globalFitInds
        [paramsEXDGlob{ss, cc}, stdEXDGlob{ss, cc}, ...
            c95EXDGlob{ss, cc}] = densityMax(rndEXDglobal{ss, cc});
        [paramsDDMGlob{ss, cc}, stdDDMGlob{ss, cc}, ...
            c95DDMGlob{ss, cc}] = densityMax(rndDDMglobal{ss, cc});

        pIndsTemp = paramsIndsNoBIC{1+cc};
        if ~loadedFlag(ss, 1)
            switch cc
                case 1
                    pIndsTemp([1, 2]) = [];
                case 2
                    pIndsTemp(2) = [];
                case 3
                    pIndsTemp(1) = [];
            end
        end
                    
        fid = fopen([saveDIR, 'DDM ', strrep(globalTitles{1+cc}, ':', ''), ...
            ' ', sources{ss}, '.txt'], 'w');
        fprintf(fid, '%s\t', paramsTitles{pIndsTemp});
        fprintf(fid, '\n');
        fprintf(fid, [repmat('%.4e\t', 1, size(rndDDMglobal{ss, cc}, 2)), '\n'], ...
            rndDDMglobal{ss, cc}');
        fclose(fid);
    end
end

nFits = nGlobalFitOpts+1;
outputString = cell(nFits, 2);
outputString(:, 1) = {sprintf('%s\n', 'Extrema detection')};
outputString(:, 2) = {sprintf('%s\n', 'DDM')};

nTabs = [...
    1, 4; ...
    1, 2; ...
    2, 1; ...
    1, 1; ...
    1, 0];

errDesc = {'std', '2.5 percentile', '97.5 percentile'};

for EXDDDM = 1:2
    if EXDDDM == 1
        numsIndep = {paramsEXDindep, stdEXDindep, ...
            cellfun(@(x) x(2:end, :), c95EXDindep, 'UniformOutput', false), ...
            cellfun(@(x) x(1:end-1, :), c95EXDindep, 'UniformOutput', false)};
        numsGlob = {paramsEXDGlob, stdEXDGlob, ...
            cellfun(@(x) x(2:end, :), c95EXDGlob, 'UniformOutput', false), ...
            cellfun(@(x) x(1:end-1, :), c95EXDGlob, 'UniformOutput', false)};
        BICInd = BICEXDIndep;
        BICGlob = BICEXDGlob;
    else
        numsIndep = {paramsDDMindep, stdDDMindep, ...
            cellfun(@(x) x(2:end, :), c95DDMindep, 'UniformOutput', false), ...
            cellfun(@(x) x(1:end-1, :), c95DDMindep, 'UniformOutput', false)};
        numsGlob = {paramsDDMGlob, stdDDMGlob, ...
            cellfun(@(x) x(2:end, :), c95DDMGlob, 'UniformOutput', false), ...
            cellfun(@(x) x(1:end-1, :), c95DDMGlob, 'UniformOutput', false)};
        BICInd = BICDDMIndep;
        BICGlob = BICDDMGlob;
    end
        
for ff = 1:nFits
    tempTabs1 = repmat(char(9), 1, nTabs(ff, 1));
    tempTabs2 = repmat(char(9), 1, nTabs(ff, 2));
    
    temp = [intFolders; repmat({tempTabs2}, 1, 3)];
    outputString{ff, EXDDDM} = [outputString{ff, EXDDDM}, ...
        sprintf('%s%s%s%s%s%s%s%s\n', globalTitles{ff}, tempTabs1, temp{:})];
    outputString{ff, EXDDDM} = [outputString{ff, EXDDDM}, ...
        sprintf('\t%s', paramsTitles{paramsIndsMaps{ff}})];
    outputString{ff, EXDDDM} = [outputString{ff, EXDDDM}, newline()];
    for ss = GTRange
        for dd = 1:numel(numsIndep)
            paramsIndep = numsIndep{dd};
            paramsGlob = numsGlob{dd};
            
            if dd == 1
                rowTitle = sources{ss};
            else
                rowTitle = errDesc{dd-1};
            end
            outputString{ff, EXDDDM} = [outputString{ff, EXDDDM}, ...
                sprintf('%s\t', rowTitle)];
            
            switch ff
                case 1
                    for ii = intRange
                        if loadedFlag(ss, ii)
                            outputString{ff, EXDDDM} = [outputString{ff, EXDDDM}, ...
                                sprintf('%.4e\t', paramsIndep{ss, ii})];
                            if dd == 1
                                outputString{ff, EXDDDM} = [outputString{ff, EXDDDM}, ...
                                    sprintf('%.0f\t', BICInd{ss}(ii))];
                            else
                                outputString{ff, EXDDDM} = [outputString{ff, EXDDDM}, ...
                                    sprintf('\t')];
                            end
                        else
                            outputString{ff, EXDDDM} = [outputString{ff, EXDDDM}, ...
                                sprintf('\t%s', tempTabs2)];
                        end
                    end
                case 2
                    if loadedFlag(ss, 1)
                        outputString{ff, EXDDDM} = [outputString{ff, EXDDDM}, ...
                            sprintf('%.4e\t', paramsGlob{ss, ff-1})];
                    else
                        outputString{ff, EXDDDM} = [outputString{ff, EXDDDM}, tempTabs2, ...
                            sprintf('%.4e\t', paramsGlob{ss, ff-1})];
                    end
                case 3
                    if loadedFlag(ss, 1)
                        outputString{ff, EXDDDM} = [outputString{ff, EXDDDM}, ...
                            sprintf('%.4e\t', paramsGlob{ss, ff-1})];
                    else
                        outputString{ff, EXDDDM} = [outputString{ff, EXDDDM}, ...
                            sprintf('%.4e\t\t', paramsGlob{ss, ff-1}(1)), ...
                            sprintf('%.4e\t', paramsGlob{ss, ff-1}(2:end))];
                    end
                case 4
                    if loadedFlag(ss, 1)
                        outputString{ff, EXDDDM} = [outputString{ff, EXDDDM}, ...
                            sprintf('%.4e\t', paramsGlob{ss, ff-1})];
                    else
                        outputString{ff, EXDDDM} = [outputString{ff, EXDDDM}, tempTabs2, ...
                            sprintf('%.4e\t', paramsGlob{ss, ff-1})];
                    end
                case 5
                    outputString{ff, EXDDDM} = [outputString{ff, EXDDDM}, ...
                        sprintf('%.4e\t', paramsGlob{ss, ff-1})];
            end
            if ff > 1 && dd == 1
                outputString{ff, EXDDDM} = [outputString{ff, EXDDDM}, ...
                    sprintf('%.0f\t', BICGlob(ss, ff-1))];
            end
            outputString{ff, EXDDDM} = [outputString{ff, EXDDDM}, newline()];
        end
        outputString{ff, EXDDDM} = [outputString{ff, EXDDDM}, newline()];
    end
end
end

fid = fopen([saveDIR, 'Extrema ', 'peaks and errors spreadsheet.txt'], 'w');
fprintf(fid, '%s\n\n', outputString{:, 1});
fclose(fid);

fid = fopen([saveDIR, 'DDM ', 'peaks and errors spreadsheet.txt'], 'w');
fprintf(fid, '%s\n\n', outputString{:, 2});
fclose(fid);

function [params, stdOut, c95] = densityMax(samples)
% To find the peak of a discretised pdf. Taking the expectation value is
% not sufficient because it does not account for skewed, or bimodal
% distributions. To find the peak, conceptually, each sample point is
% convolved with a Gaussian before summing together and finding the maximum 
% of the smoothed pdf. 
%
% For general pdfs, choosing too small a width of the Gaussian will bias
% the result toward small clusters of closely spaced samples. A very large
% width beomes equivalent to taking the expectation value. The 
% multi-dimensional Gaussian width is the range of the diagonal multiplied
% by the rangeFactor.
% 
% However, for the data in the decision making publication (Wong et al), 
% the peak calculation is quite insensitive to the rangeFactor. I tried a
% range of values from 1/100 to 100, with negligible changes in the output.
% Copyright (C) 2023 Clifford Talbot

rangeFactor = 3; % this is the value used for the analysis presented.
[nData, nParams] = size(samples);

% compute range of data and Gaussian widths
UU = max(samples, [], 1);
LL = min(samples, [], 1);
init = LL + (UU-LL)/2;
sigmas = (UU - LL)*nParams*rangeFactor; % Gaussian width = sigmas^2

% find the maximum of the smoothed pdf using curve minimisation algorithm
prob([], sigmas, samples); % initialise helper function
hFn = @(x) 1/prob(x, [], []); % want to maximise prob = minimise 1/prob
opts = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', ...
    'FunctionTolerance', 1e-12, 'MaxFunctionEvaluations', 600*nParams, ...
    'MaxIterations', 1000, 'StepTolerance', 1e-12, 'Display', 'off');
[params, ~, ~, exitFlag] = lsqnonlin(hFn, init, [], [], opts); 
if exitFlag < 0
    disp('Error finding max density');
end

% compute spread in the non-smoothed pdf (standard deviation and 95 % CI)
centred = samples - ones(size(samples, 1), 1)*params;
stdOut = sqrt(mean(centred.*centred));
c95 = zeros(2, nParams);
uuFr = nData*0.975;
uuT1 = floor(uuFr); uuT2 = ceil(uuFr); uuT3 = uuFr-uuT1;
hUU = @(x) ((x(uuT2)-x(uuT1))*uuT3) + x(uuT1);
llFr = nData*0.275;
llT1 = floor(llFr); llT2 = ceil(llFr); llT3 = llFr-llT1;
hLL = @(x) ((x(llT2)-x(llT1))*llT3) + x(llT1);
for pp = 1:nParams
    sorted = sort(samples(:, pp));
    c95(1, pp) = hUU(sorted);
    c95(2, pp) = hLL(sorted);
end
end

function p = prob(coords, sigmasIn, samplesIn)
% compute the value of the smoothed pdf at point defined by coords
persistent sigmas samples nData tempOnes sigmaSq const

if isempty(coords)
    sigmas = sigmasIn;
    samples = samplesIn;
    [nData, nCoords] = size(samples);
    tempOnes = ones(nData, 1);
    sigmaSq = 1./(2 * tempOnes .* (sigmas.*sigmas));
    const = (1/nData) * prod(1./sigmas) * ((1/(2*pi))^(nCoords/2));
    return;
end

shiftSamples = samples - tempOnes*coords;
scaledDelta = sum(...
    shiftSamples.*shiftSamples .* sigmaSq, 2);
p = const * sum(exp(-scaledDelta));
end
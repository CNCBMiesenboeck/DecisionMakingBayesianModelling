function p = pDDM(reactionTimes, choiceAccuracy, difficulty, A, K, tND, sigmaND)
% to compute the log-probability of a reaction time and choice accuracy at
% a task difficulty level given the model parameters A, k, tND, sigmaND for
% the drift diffusion model. 
% Copyright (C) 2023 Clifford Talbot

if ~checkBounds(A, K, tND, sigmaND)
    p = -inf;
    return;
end

err = 1e-29; % numerical error limit for Navarro-Fuss approximation

% solution to DDM has bounds at 0 and 1. Convert our A and K appropriately
a = 2*A;
v = 2*K * difficulty;

% the distribution of the non-decision time is numerically convolved with
% the drift diffusion model. Each measured reation time is effectively a
% delta function, so the convolution becomes a sumation of Gaussian 
% functions. Each reactionTime is replaced with a series of time points 
% sampling the range of the Gaussian (Gaussian-sampled reaction times). The
% DDM probability for each of these time points is computed and weighted by
% the Gaussian function. The log-probabilities are then summed.
nTSteps = 31; % number of points in the non-decision time distribution
tMaxSD = 5*sigmaND; % use a time range of 5 standard deviations
dt = 2*tMaxSD/(nTSteps-1);
tNDRange = ((-tMaxSD):dt:tMaxSD); % the timepoints sampling the Gaussian (centred at zero)

nData = numel(reactionTimes);
nData_RTND = nData*nTSteps; % total number of Gaussian-sampled reaction times (RT=reactionTime, ND=non-decision)
onesNT = ones(1, nTSteps); % convenient arrays
onesDD = ones(nData, 1);

c1 = dt/(sigmaND*sqrt(2*pi)); % Gaussian coefficients
c2 = -1/(2*sigmaND*sigmaND);
hG = @(x) c1*exp(x.*x * c2); % Gaussian function
sigmaNDWeighting = hG(tNDRange); 
weightingMat = onesDD*sigmaNDWeighting; % each row is a zero-centred Gaussian

t_RTNDraw = reactionTimes*onesNT + (onesDD*tNDRange) - tND; % each row 
% contains the Gaussian-sampled reaction times for each measured reaction
% time

% compute the minimum Gaussian-sampled reaction times for each measured
% reactionTime. The Gaussian-sampled time range is shifted forwards for any
% non-physical negative values.
tMins = (reactionTimes-(tND+tMaxSD));
tSteps = (1:nTSteps) * dt;
for ii = find(tMins<=0)'
    t_RTNDraw(ii, :) = tSteps;
    weightingMat(ii, :) = hG(tNDRange-tMins(ii));
end

% reshape data, and replicate if necessary, for passing to wfpt
t_RTND = reshape(t_RTNDraw, nData_RTND, 1);
weighting = reshape(weightingMat, nData_RTND, 1);
ch_RTND = reshape(choiceAccuracy*onesNT, nData_RTND, 1);
v_RTND = reshape(v*onesNT, nData_RTND, 1);

p_RTND = zeros(nData_RTND, 1); % initialise output of wfpt
tCheck = t_RTND>0;
chT = ch_RTND & tCheck; % correct choices
chBar = (~ch_RTND) & tCheck; % incorrect
p_RTND(chT) = wfpt(t_RTND(chT), -v_RTND(chT), a, err, false); % compute probabilities of correct choices
p_RTND(chBar) = wfpt(t_RTND(chBar), v_RTND(chBar), a, err, false); % incorrect choices
p = log(sum(reshape(p_RTND.*weighting, nData, nTSteps), 2)); % apply Gaussian weighting and sum over the Gaussian
end
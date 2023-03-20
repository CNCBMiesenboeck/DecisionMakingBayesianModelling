function p = logP_RT_acc(reactionTimes, correctChoice, difficulty, A, K, tND, sigmaND, dt)
% to compute the log-probability of a reaction time and choice accuracy at
% a task difficulty level given the model parameters A, k, tND, sigmaND and
% dt for the extrema detection model. The equations are as defined in the
% methods section of the publication.
% Copyright (C) 2023 Clifford Talbot

if ~checkBounds(A, K, tND, sigmaND)
    p = -Inf;
    return;
end

temp1 = sqrt(2*dt);
pPos = 0.5*erfc((A-K*difficulty*dt)/temp1);
pNeg = 0.5*erfc((A+K*difficulty*dt)/temp1);
pTot = pPos + pNeg;

pNormPos = pPos./pTot;

t1 = pTot/dt;
t2 = tND - reactionTimes;
t3 = sigmaND*t1;
sr2 = sqrt(2);
t4 = t3/sr2;
t5 = 1/(sigmaND*sr2);
t6 = t2*t5 + t4;

p = log(0.5*t1) + ((t2 + t3*sigmaND/2) .* t1);
p = p + log(erfc(t6));
p(correctChoice) = p(correctChoice) + log(pNormPos(correctChoice));
chBar = ~correctChoice;
p(chBar) = p(chBar) + log(1-pNormPos(chBar));
end
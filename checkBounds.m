function TF = checkBounds(A, K, tND, sigmaND)
% return true if in bounds
% It is not physical for any parameter to be negative
% It is not realistic for the non-decision time, or it's spread to be 
% greater than 3 seconds
TF = ~((A <= 0) || (K <= 0) || (tND <= 0) || (tND > 3) || ...
        (sigmaND > 3) || ...
        (sigmaND < 0)); 
end
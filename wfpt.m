function p=wfpt(t, v, a, err, logFlag)
% code from D. J. Navarro and I. G. Fuss. Journal of Mathematical 
% Psychology 53 (2009) 222-230. With minor modifications to:
%   accept multiple data points more efficiently than a simple loop
%   compute log-probability
%   always start from midpoint between upper and lower bounds

nTime = numel(t);
p = zeros(nTime, 1); % initialise output array

tt = t/(a^2); % use normalized time
w = 0.5; % start at midpoint between bounds (changed from w=z/a in Navarro&Fuss)
for ii = 1:nTime
    p(ii) = navarroCode(tt(ii), w, err);
end
% convert to f(t|v,a,w)
if logFlag
    p = log(p) + (-v*a*w -(v.^2).*t/2) - log(a^2);
else
    p=p.*exp(-v*a*w -(v.^2).*t/2)/(a^2);
end
end

function p = navarroCode(tt, w, err)
% directly from Navarro and Fuss

% calculate number of terms needed for large t
if pi*tt*err<1 % if error threshold is set low enough
    kl=sqrt(-2*log(pi*tt*err)./(pi^2*tt)); % bound
    kl=max(kl,1/(pi*sqrt(tt))); % ensure boundary conditions met
else % if error threshold set too high
    kl=1/(pi*sqrt(tt)); % set to boundary condition
end

% calculate number of terms needed for small t
if 2*sqrt(2*pi*tt)*err<1 % if error threshold is set low enough
    ks=2+sqrt(-2*tt.*log(2*sqrt(2*pi*tt)*err)); % bound
    ks=max(ks,sqrt(tt)+1); % ensure boundary conditions are met
else % if error threshold was set too high
    ks=2; % minimal kappa for that case
end
% compute f(tt|0,1,w)
p=0; %initialize density
if ks<kl % if small t is better...
    K=ceil(ks); % round to smallest integer meeting error
    for k=-floor((K-1)/2):ceil((K-1)/2) % loop over k
        p=p+(w+2*k)*exp(-((w+2*k)^2)/2/tt); % increment sum
    end
    p=p/sqrt(2*pi*tt^3); % add constant term
else % if large t is better...
    K=ceil(kl); % round to smallest integer meeting error
    for k=1:K
        p=p+k*exp(-(k^2)*(pi^2)*tt/2)*sin(k*pi*w); % increment sum
    end
    p=p*pi; % add constant term
end
% % convert to f(t|v,a,w)             Do this outside the loop (avoids 
% p=p*exp(-v*a*w -(v^2)*t/2)/(a^2);   repeatedly computing the second term)
end
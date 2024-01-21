%% Section 2.1. Preprocessing (median_filter)
%%
% In order to reduce the effect of ectopic beats in the RR series, we may
% use a simple 3-point median filter as: rm(n)=median{r(n-1),r(n),r(n+1)}.
% This filter is also useful to reject the outlier RR intervals.

function rm = median_filter(r)
    rm = medfilt1(r);
end

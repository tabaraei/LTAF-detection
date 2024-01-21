%% Section 2.3. Bigeminy Supression (bigeminy_supressor)
%%
% Bigeminy is a cardiac arrhythmia characterized by a pattern in which 
% every normal heartbeat is followed by a premature beat, creating a 
% regular pattern of two beats, which can be incorrectly interpreted as AF
% when the detection is RR-based.

function [B, Bt] = bigeminy_supressor(r, rm, N_intervals, N, alpha)
    
    % Computing B(n)
    B = zeros(size(r));
    for n=1:N_intervals
        start_idx = max(1, n - (N-1));
        numerator = sum(rm(start_idx:n));
        denominator = sum(r(start_idx:n));
        B(n) = ((numerator/denominator) - 1) ^ 2;
    end

    % Computing the smoothed version of B(n)
    Bt = forward_backward_averager(B, alpha);

end


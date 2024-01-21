%% Section 2.2. RR Irregularity Detection (irregularity_detector)
%%
% In order to distinguish RR irregularities, we should use a sliding 
% detection window of length N, located at time n, and compute the number 
% of all pairwise RR interval combinations differing more than gamma 
% seconds, and normalize them with their maximum value. To achieve this:
%
% 1. First, we will introduce the heaviside function "H", where it simply
% denotes whether the pairwise difference of given RR intervals are below a
% certain threshold "gamma" or not.
%
% 2. For each RR interval "n", we will compute the count of pairwise RR
% interval combinations differing more than "gamma", and we will normalize
% this value w.r.t its maximum value.
%
% 3. We compute the ratio "I" between the smoothed version of M utilizing
% exponential averaging, and the RR interval trend stored in "rt".

function [M, Mt, It] = irregularity_detector(rm, rt, N_intervals, N, gamma, alpha)

    % Defining function H(.) as a Heaviside Step function
    H = @(r1, r2) double(abs(r1 - r2) >= gamma);
    
    % Computing the M(n) as a sliding window of length N
    M = zeros(size(rm));
    for n=1:N_intervals
        count = 0;
        for j=0:N-2
            for k=j+1:N-1
                if (n-j)>0 && (n-k)>0
                    count = count + H(rm(n-j), rm(n-k));
                end
            end
        end
        M(n) = count * 2/(N*(N-1));
    end
    
    % Computing the smoothed version of M(n)
    Mt = forward_backward_averager(M, alpha);
    
    % Computing the RR irregularity ratio denoted as I
    It = Mt ./ rt;

end
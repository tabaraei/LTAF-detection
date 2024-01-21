%% Section 2.1. Preprocessing (forward_backward_averager)
%%
% This function computes the forward-backward filtering to achieve a linear
% phase on a given series X, and it aims to compute the exponential average
% to better track the trend using 0<alpha<1 as the degree of smoothing.
% This procedure results in the estimation of the mean RR interval, which 
% can be used as a feature in the AF detector.
%
% We use the filtfilt Zero-phase forward and reverse digital IIR filtering.
%
%     Y = filtfilt(B, A, X) filters the data in vector, matrix, or N-D
%     array, X, with the filter described by vectors A and B to create
%     the filtered data Y. The filter is described by the difference
%     equation:
%  
%       a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                             - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
%     
%       
%     In order to achieve the desired behaviour by the paper, which is:
% 
%       y(n) = (1-alpha)*y(n-1) + alpha*x(n)
%     
%     We have to set the following coefficients in the difference equation:
%
%       a(1) = 1
%       a(2) = alpha - 1
%       b(1) = alpha

function Y = forward_backward_averager(X, alpha)
    A = [1, alpha-1];
    B = alpha;
    Y = filtfilt(B, A, X);
end

%% Section 2.4. Signal Fusion and Detection (signal_fusion)
%%
% The decision function O(n) is produced through a simple signal fusion,
% which is identical to "Bt" unless it exceeds a fixed threshold delta
% when, instead it becomes identical to "It".

function [O, RR_AF_Predictions] = signal_fusion(It, Bt, delta, eta)

    % Calculation of the O(n)
    O = It.*(Bt >= delta) + Bt.*(Bt < delta);

    % Finding the AF where O(n) is above the eta threshold
    RR_AF_Predictions = O > eta;

end


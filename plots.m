%% Plotting Functions
%%
% This file contains written codes to plot a specific diagram, and they are
% organized for a more clean code

function [plt] = plots()
    
    % function handlers accessible in the "main.m" file
    plt = struct();
    plt.plot_whole_signal = @pws;
    plt.plot_signal_channels = @psc;
    plt.plot_rr_intervals = @prr;
    plt.plot_median_filtering = @pmf;
    plt.plot_exponential_averager = @pea;
    plt.plot_irregularity_detector = @pid;
    plt.plot_bigeminy_supressor = @pbs;
    plt.plot_signal_fusion = @psf;
    plt.plot_final_results = @pfr;

    % plotting colors
    color_primary = "#4494bd";
    color_secondary = "#bd444c";
    color_tertiary = "#55a868";

    %% pa - Plotting the annotations
    %%
    % This function takes as an input a vector x as RR series, and plot the
    % annotations on top of it. This annotations are however not so useful,
    % so the function is no further used in the algorithm.
    %
    function pa(x, annots, text)
        unique_annots = unique(annots);
        legends = cell(1, length(unique_annots));
        annotationMapping = containers.Map(...
            {'N', '|', 'T', 'V', 'A', 'Q', '+', '"'}, ...
            {'Normal Beat', 'Detected Artifacts', 'AF Termination', ...
            'Premature Ventricular Contraction', ...
             'Atrial Premature Beat', 'Unclassifiable Beat', ...
             'Rhythm Change', 'Comment Annotation'});

        hold on;
        plot(x);
        legends{1} = text;
        idx = 2;
        for i = 1:length(unique_annots)
            annot_indices = find(annots == unique_annots(i));
            if annot_indices ~= 'N'
                scatter(annot_indices, x(annot_indices), 20, 'filled');
                legends{idx} = annotationMapping(unique_annots(i));
                idx = idx + 1;
            end
        end
        hold off;
        title([text, ' and their Annotations']);
        legend(legends);
        xlabel('Index');
        ylabel('RR Interval');
    end

    %% pax - Plotting the auxiliary annotations
    %%
    % Similar to the function "pa" above, which instead uses the auxiliary
    % annotations (comments), where it contains the starting points of each
    % episode in the signal
    function pax(x, annots_aux, text)
        unique_annots = unique(annots_aux(~cellfun('isempty', annots_aux)));
        legends = cell(1, length(unique_annots));
        annotationMapping = containers.Map(...
            {'(N', '(SVTA', '(VT', '(AFIB', '(B', '(T', '(IVR', '(AB', ...
            '(SBR', ' Aux', 'MISSB', 'PSE', 'M', 'MB'}, ...
            {'Normal sinus rhythm', 'Supraventricular tachyarrhythmia', ...
            'Ventricular tachycardia', 'Atrial fibrillation', ...
            'Ventricular bigeminy', 'Ventricular trigeminy', ...
            'Idioventricular rhythm', 'Atrial bigeminy', ...
            'Sinus bradycardia', 'Auxiliary', 'Missed beat', ...
            'Pause', 'M', 'MB'});

        hold on;
        plot(x, 'Color', color_primary);
        legends{1} = text;
        for i = 1:length(unique_annots)
            annot_indices = find(strcmp(annots_aux, unique_annots{i}));
            scatter(annot_indices, x(annot_indices), 20, 'filled');
            legends{i+1} = annotationMapping(unique_annots{i});
        end
        hold off;
        title([text, ' and Their Starting Episode Annotations']);
        legend(legends);
        ylabel('R-R(t) in seconds');
    end


    %% pws - Plotting the whole signal
    %%
    % Plotting the whole ECG signal with all its channels, while marking 
    % the AF episodes on the plot
    %
    function pws(signal, ECG_AF_Groundtruth, N_channels)
        figure()
        for channel=1:N_channels
            % create the Atrial Fibrilation mask
            AF_mask = signal(:,channel).*ECG_AF_Groundtruth;
            AF_mask_indices = find(AF_mask~=0);

            % Plot ECG channels and mask with AF
            subplot(N_channels, 1, channel);
            hold on;
            plot(signal(:,channel), 'Color', color_primary);
            plot(AF_mask_indices, AF_mask(AF_mask_indices), 'Color', color_secondary)
            hold off;
            ylabel('Amplitude (mV)');
            if sum(ECG_AF_Groundtruth) > 0
                legend('ECG Signal', 'Atrial Fibrilation Episodes');
            else
                legend('ECG Signal');
            end
            
            title(['Channel ', num2str(channel)]);
        end
        sgtitle('Complete ECG Signal');
    end

    %% psc - Plotting the signal channels and annotation indices
    %%
    % Plotting the original signal with all its channels, in a 10-second 
    % period, and add the corresponding annotations of QRS indices where
    % they exist (middle point)
    %
    function psc(signal, indices, N_channels, N_intervals, fs)
        start_plot = indices(floor(N_intervals/2)) - fs;
        end_plot = start_plot + 10*fs;

        figure()
        sgtitle('Different ECG Channels and Position of Annotations');
        for channel=1:N_channels
            subplot(N_channels, 1, channel);
            hold on;
            plot(signal(:,channel), 'Color', color_primary);
            plot(indices, signal(indices, channel), '*', 'Color', color_secondary);
            hold off;
            ylabel('Amplitude (mV)');
            xlim([start_plot, end_plot]);
            legend(['Channel ', num2str(channel)], 'QRS Annotation');
            title(['Channel ', num2str(channel)]);
        end
    end

    %% prr - Plotting the RR intervals
    %%
    function prr(r, annots_aux)
        figure()
        sgtitle('RR Intervals and their Annotations')

        subplot(2,1,1);
        plot(r, 'Color', color_primary);
        title('RR Intervals');
        ylabel('R-R(t) in seconds');

        subplot(2,1,2);
        pax(r, annots_aux, 'RR Intervals')
        
    end

    %% pmf - Plotting the 3-point median filtering results
    %%
    function pmf(r, rm)
        figure()
        hold on;
        plot(r, 'Color', color_primary);
        plot(rm, 'Color', color_secondary);
        hold off;
        legend('Original RR Series', 'Filtered RR Series');
        ylabel('R-R(t) in seconds');
        title('Filtering RR Intervals with 3-Point Median Filter');
    end

    %% pea - Plotting the forward-backward exponential averager
    %%
    function pea(r, rt)
        figure()
        hold on;
        plot(r, 'Color', color_primary);
        plot(rt, 'LineWidth', 2, 'Color', color_secondary);
        hold off;
        legend('Filtered RR Series', 'Averaged RR Series');
        ylabel('R-R(t) in seconds');
        title('Forward-Backward Exponential Averaged RR Intervals');
    end

    %% pid - Plotting the irregularity detection
    %%
    function pid(M, Mt, It)
        figure()

        subplot(2,1,1);
        hold on;
        plot(M, 'Color', color_primary);
        plot(Mt, 'LineWidth', 2, 'Color', color_secondary);
        hold off;
        ylabel('M(n)');
        legend('M(n)', 'Smoothed M(n)');
        title('M(n) vs. Forward-Backward Exponential Average of M(n)');

        subplot(2,1,2);
        plot(It, 'Color', color_primary);
        ylabel('I(n)');
        title('Primary Feature of RR Irregularity Denoted as "I"');
        annotation('textbox', [0.15, 0.25, 0.18, 0.18], ...
            'String', 'I(n) having higher values during AF');
        annotation('textbox', [0.7, 0.15, 0.18, 0.18], ...
            'String', 'I(n) close to 0 for regular rhythms');
    end

    %% pbs - Plotting the bigeminy supression
    %%

    function pbs(B, Bt, N_intervals, delta)
        figure()
        hold on;
        plot(B, 'Color', color_primary);
        plot(Bt, 'LineWidth', 2, 'Color', color_secondary);
        line([1, N_intervals], [delta, delta], 'Color', 'green', 'LineStyle', '--');
        hold off;
        ylabel('B(n)');
        legend('B(n)', 'Smoothed B(n)', 'Delta Threshold');
        title('Bigeminy Suppression Measure of RR Irregularity Denoted as "B"');
        annotation('textbox', [0.2, 0.6, 0.2, 0.2], ...
            'String', 'B(n) might indicate AF irregularity in variations');
        annotation('textbox', [0.6, 0.2, 0.2, 0.2], ...
            'String', 'B(n) close to 0 for regular or bigeminy rhythms');
    end

    %% psf - Plotting the signal fusion and detection
    %%
    function psf(O, annots_aux, N_intervals, RR_AF_Predictions, RR_AF_Groundtruth, eta)
        figure()

        % O(n) vs. predicted atrial fibrillation
        AF_preds_mask = O.*RR_AF_Predictions;
        AF_preds_mask_indices = find(AF_preds_mask~=0);
        subplot(3,1,1);
        hold on;
        plot(O, 'Color', color_primary);
        if sum(RR_AF_Predictions) > 0
            plot(AF_preds_mask_indices, AF_preds_mask(AF_preds_mask_indices), 'Color', color_secondary);
            line([1, N_intervals], [eta, eta], 'Color', color_tertiary, 'LineStyle', '--');
            legend('O(n)', 'Predicted Atrial Fibrillation', 'Eta Threshold');
        else
            line([1, N_intervals], [eta, eta], 'Color', color_tertiary, 'LineStyle', '--');
            legend('O(n)', 'Eta Threshold');
        end
        hold off;
        ylabel('O(n)');
        title('Decision Function O(n) and Predicted AF');
        ylim_copy = ylim;

        % O(n) vs. groundtruth atrial fibrillation
        AF_groundtruth_mask = O.*RR_AF_Groundtruth;
        AF_groundtruth_mask_indices = find(AF_groundtruth_mask~=0);
        subplot(3,1,2)
        hold on;
        plot(O, 'Color', color_primary);
        if sum(RR_AF_Groundtruth) > 0
            plot(AF_groundtruth_mask_indices, AF_groundtruth_mask(AF_groundtruth_mask_indices), 'Color', color_secondary);
            line([1, N_intervals], [eta, eta], 'Color', color_tertiary, 'LineStyle', '--');
            legend('O(n)', 'Groundtruth Atrial Fibrillation', 'Eta Threshold');
        else
            line([1, N_intervals], [eta, eta], 'Color', color_tertiary, 'LineStyle', '--');
            legend('O(n)', 'Eta Threshold');
        end
        hold off;
        ylabel('O(n)');
        title('Decision Function O(n) and Groundtruth AF');

        % O(n) vs. starting episode annotations
        subplot(3,1,3);
        pax(O, annots_aux, 'O(n)');
        ylim(ylim_copy);
    end

    %% pfr - Plotting the final results for first ECG channel
    %%
    % Plotting the resulting predicted AF vs. groundtruth AF for a given
    % channel.
    %
    function pfr(signal, ECG_AF_Groundtruth, indices, RR_AF_Predictions)
        channel = 1;

        % plot the ECG groundtruth mask
        AF_groundtruth_mask = signal(:,channel).*ECG_AF_Groundtruth;
        AF_groundtruth_mask_indices = find(AF_groundtruth_mask~=0);

        % plot the ECG predicted AF mask
        AF_pred_indices = indices(find(RR_AF_Predictions~=0));
        AF_pred_mask = signal(AF_pred_indices,channel);

        figure()
        hold on;
        plot(signal(:,channel), 'Color', color_primary);
        plot(AF_groundtruth_mask_indices, AF_groundtruth_mask(AF_groundtruth_mask_indices), 'Color', color_secondary)
        plot(AF_pred_indices, AF_pred_mask, '.', 'MarkerSize', 5, 'Color', color_tertiary)
        hold off;
        ylabel('Amplitude (mV)');
        legend('ECG Signal', 'Groundtruth AF Episodes', 'Predicted AF Episodes');
        title('Predicted vs. Groundtruth Atrial Fibrilation of the ECG');

    end

end
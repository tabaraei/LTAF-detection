clc, clear, close all;
addpath(fullfile(pwd, 'packages/WFDB'));
addpath(fullfile(pwd, 'functions'));
savepath;
[plt] = plots();
[LTAF] = load_data();

%% ALGORITHM
%%
% This .m file loads a single piece of data with its annotations, and runs
% the whole algorithm proposed by the paper and provides visualizations. We
% can use either "atr" or "qrs" annotation as the <annot_type>.
%
% There are some parameters investigated in the results section of the
% paper, which we will define as below:
%
%  alpha: smoothing factor of the exponential average (0<alpha<1)
%  N: size of the sliding window to detect pairwise RR interval differences
%  gamma: threshold in seconds for detection of pairwise differences in M
%  delta: threshold employed in the signal fusion and detection
%  eta: Final decision threshold
%
% Using the parameters below, we can control the algorithm behaviour when
% running the algorithm, to display a certain information, and whether to
% only consider a single file (for test purposes) or run the algorithm on
% the whole dataset:
%
% disp_annot_counts: if "true", displays the unique annotations in console
% disp_plots: if "true", displays various plots for the domain
% disp_evaluation: if "true", calculates evaluation metrics on the dataset
% test_single_file: if "true", only runs the algorithm on a single file 
% specified in the "test_file" variable

% Hyperparameters of the paper
alpha = 0.02;
N = 8;
gamma = 0.03;
delta = 2e-4;
eta = 0.725;

% Data directory and the annotation type
data_directory = 'data/LTAF/';
annot_type = 'atr';
test_file = 'data/LTAF/00';
data_paths = LTAF.get_data_paths(data_directory, annot_type);

% display settings
test_single_file = false;
disp_annot_counts = true;
disp_plots = true;
disp_evaluation = true;

% evaluation metrics
cols = {'Record', 'Accuracy', 'Precision', 'Sensitivity', 'Specificity'...
    'ECG_time', 'AF_time', 'N_AF_Episodes', 'Runtime'};
evaluation_metrics = array2table(...
    strings(length(data_paths), length(cols)), 'VariableNames', cols);

% start of the algorithm
for i=1:length(data_paths)

    % set the path as "test_file" if we don't want to loop over all files
    if test_single_file
        path = test_file;
    else
        path = data_paths{i};
    end
    disp(repmat('-', 1, 80));
    disp(['File name: ', path]);

    % Get the data and groundtruth AF annotations
    [signal, indices, annots, annots_aux, r, N_channels, N_intervals, fs] = ...
        LTAF.get_data(path, annot_type, disp_annot_counts);
    [ECG_AF_Groundtruth, RR_AF_Groundtruth, ECG_time, AF_time, N_AF_Episodes] = ...
        LTAF.get_annots(signal, indices, annots_aux, N_intervals, fs);
    disp(['Whole ECG duration: ', ECG_time]);
    disp(['AF episodes duration: ', AF_time]);
    
    % Main algorithm, computing the runtime within (tic, toc)
    tic;
    rm = median_filter(r);
    rt = forward_backward_averager(r, alpha);
    [M, Mt, It] = irregularity_detector(rm, rt, N_intervals, N, gamma, alpha);
    [B, Bt] = bigeminy_supressor(r, rm, N_intervals, N, alpha);
    [O, RR_AF_Predictions] = signal_fusion(It, Bt, delta, eta);
    runtime = toc;

    % display the according plots
    if disp_plots

        % Plotting the signal channels and RR intervals
        plt.plot_whole_signal(signal, ECG_AF_Groundtruth, N_channels)
        plt.plot_signal_channels(signal, indices, N_channels, N_intervals, fs)
        plt.plot_rr_intervals(r, annots_aux)

        % 3-Point Median Filtering
        plt.plot_median_filtering(r, rm)

        % Forward-Backward Exponential Averaging
        plt.plot_exponential_averager(r, rt)
        
        % RR Irregularity Detection
        plt.plot_irregularity_detector(M, Mt, It)
        
        % Bigeminy Supression
        plt.plot_bigeminy_supressor(B, Bt, N_intervals, delta)

        % Signal Fusion and Detection
        plt.plot_signal_fusion(O, annots_aux, N_intervals, RR_AF_Predictions, RR_AF_Groundtruth, eta)
        plt.plot_final_results(signal, ECG_AF_Groundtruth, indices, RR_AF_Predictions)

    end

    % store the performance evaluations
    if disp_evaluation
        % computing confusion matrix and extracting values
        % confMat = confusionmat(RR_AF_Groundtruth, RR_AF_Predictions);
        TN = sum(~RR_AF_Groundtruth & ~RR_AF_Predictions);
        FP = sum(~RR_AF_Groundtruth & RR_AF_Predictions);
        FN = sum(RR_AF_Groundtruth & ~RR_AF_Predictions);
        TP = sum(RR_AF_Groundtruth & RR_AF_Predictions);
        
        % evaluation metrics
        evaluation_metrics.Record{i} = path;
        evaluation_metrics.Accuracy{i} = num2str((TP + TN) / N_intervals);
        evaluation_metrics.Precision{i} = num2str(TP / (TP + FP));
        evaluation_metrics.Sensitivity{i} = num2str(TP / (TP + FN));
        evaluation_metrics.Specificity{i} = num2str(TN / (TN + FP));
        evaluation_metrics.ECG_time{i} = ECG_time;
        evaluation_metrics.AF_time{i} = AF_time;
        evaluation_metrics.N_AF_Episodes{i} = num2str(N_AF_Episodes);
        evaluation_metrics.Runtime{i} = num2str(runtime);
    end

    if test_single_file
        break;
    end
end

if disp_evaluation
    % Write the table to the CSV file
    writetable(evaluation_metrics, 'evaluations.csv');
    disp(repmat('-', 1, 80));
    disp(evaluation_metrics)
end

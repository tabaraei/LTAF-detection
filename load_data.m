%% Preparing the Data and Annotations
%%
% Long-Term Atrial Fibrillation Database Installation
% 
% WFDB Package Download URL:
%    "https://physionet.org/static/published-projects/matlab-toolbox/wfdb-app-matlab-0.10.0.zip"
%
% Data Download Link:
%    "https://physionet.org/content/ltafdb/1.0.0"
%
% Considering that the WFDB package is downloaded, unzipped, and all the
% contents of "mcode" folder is present at subfolder "packages/WFDB", we
% can use the following commands to install the package:
%    addpath('packages/WFDB');
%    savepath;
%
% QRS Annotations:
%    "https://physionet.org/content/ltafdb/1.0.0/"
%    N      All detected beats (including occasional ventricular ectopic beats)
%    |      Detected artifacts
%    T      AF terminations
%
% ATR Annotations:
%    "https://archive.physionet.org/physiobank/annotations.shtml"
%    N		Normal beat
%    V		Premature ventricular contraction
%    A		Atrial premature beat
%    Q		Unclassifiable beat
%    +		Rhythm change
%    "      Comment annotation

function [LTAF] = load_data()
    
    LTAF = struct();
    LTAF.get_data_paths = @gdp;
    LTAF.get_data = @gdata;
    LTAF.get_annots = @gannots;

    %% Get the data paths sorted by file name
    %%
    % This function will look for the annotation files of type "annot_type"
    % in the specified "data_directory" and sorts them by their file name.
    %
    function data_paths = gdp(data_directory, annot_type)
        % Sorting the files according to file_names integer value
        data_files = {dir([data_directory, '*.', annot_type]).name};
        custom_sort = @(x) sscanf(x, '%d');
        numeric_values = cellfun(custom_sort, data_files);
        [~, sorted_indices] = sort(numeric_values);
        sorted_file_names = data_files(sorted_indices);
        
        % Store the data path for each ECG file in a cell array
        data_paths = cell(1, length(sorted_file_names));
        for i=1:length(sorted_file_names)
            data_paths{i} = [data_directory, sorted_file_names{i}(1:end-4)];
        end
    end

    %% Loading the Long-Term Atrial Fibrillation Databse
    %%
    % LTAF database is locally stored in "data/LTAF", and is loaded for
    % further computation. This function also computes the RR series and
    % returns it in the variable "r".
    %
    function [signal, indices, annots, annots_aux, r, ...
            N_channels, N_intervals, fs] ...
            = gdata(file_name, annot_type, display_annot_count)

        % Loading the data
        [signal, fs] = rdsamp(file_name);
        [indices, annots,~,~,~, annots_aux] = rdann(file_name, annot_type);
        N_channels = size(signal, 2);

        % remove out-of-bound indices which exist in some files, e.g. "30"
        condition = indices <= length(signal) & indices > 0;
        indices = indices(condition);
        annots = annots(condition);
        annots_aux = annots_aux(condition);
   
        % Display the unique annotations found in the file
        if display_annot_count
            disp('Unique Annotations:');
            unique_annots = unique(annots_aux(~cellfun('isempty', annots_aux)));
            for i=1:length(unique_annots)
                count = sum(strcmp(annots_aux, unique_annots{i}));
                disp(['     Annotation "', unique_annots{i}, '" Count: ', num2str(count)]);
            end
        end
    
        % Preparing the RR-intervals same size as annotations (fs not employed)
        r = diff(indices/fs);
        r = [r(1); r];
        N_intervals = length(r);
    
    end

    %% Building the ECG and RR Interval Groundtruth Labels
    %%
    % In order to create the groundtruths for the signal and RR series, we
    % need to find the AF episodes starting from the annotation "(AFIB"
    % until the end of episode, which shows the whole AF episode and is
    % marked as "1" in the groundtruths. Wherever AF is not present, the
    % groundtruth is "0".
    %
    function [ECG_AF_Groundtruth, RR_AF_Groundtruth, ECG_time, AF_time, N_AF_Episodes] ...
            = gannots(signal, indices, annots_aux, N_intervals, fs)

        % create the begin-to-end intervals of episode-level annotations for RR series
        RR_begin_episodes = find(cellfun(@(x) ischar(x) && startsWith(x, '('), annots_aux));
        RR_end_episodes = [RR_begin_episodes(2:end)-1; N_intervals];
        
        % create the begin-to-end intervals of annotations for ECG signal
        ECG_begin_episodes = indices(RR_begin_episodes);
        ECG_end_episodes = [ECG_begin_episodes(2:end)-1; length(signal)];
        
        % create the AF groundtruths for both ECG and RR intervals
        label_episodes = annots_aux(RR_begin_episodes);
        RR_AF_Groundtruth = zeros(N_intervals, 1);
        ECG_AF_Groundtruth = zeros(length(signal), 1);
        N_AF_Episodes = 0;

        for i=1:length(label_episodes)
            if strcmp(label_episodes{i}, "(AFIB")
                N_AF_Episodes = N_AF_Episodes + 1;
                RR_AF_Groundtruth(RR_begin_episodes(i):RR_end_episodes(i)) = 1;
                ECG_AF_Groundtruth(ECG_begin_episodes(i):ECG_end_episodes(i)) = 1;
            end
        end
        ECG_AF_Groundtruth = boolean(ECG_AF_Groundtruth);
        RR_AF_Groundtruth = boolean(RR_AF_Groundtruth);
        
        % display the ECG and AF durations in custom format
        % ECG_time = datestr(seconds(length(signal) / fs), 'dd HH:MM:SS');
        % AF_time = datestr(seconds(sum(ECG_AF_Groundtruth) / fs), 'dd HH:MM:SS');
        hours = @(X) floor(X / 3600);
        minutes = @(X) floor((X - hours(X) * 3600) / 60);
        seconds = @(X) round(mod(X, 60));
        custom_time = @(X) sprintf('%02d:%02d:%02d', hours(X), minutes(X), seconds(X));
        ECG_time = custom_time(length(signal) / fs);
        AF_time = custom_time(sum(ECG_AF_Groundtruth) / fs);
    end

end

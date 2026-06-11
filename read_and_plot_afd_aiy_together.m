clear;
clc;
% close all;

%% ================= USER INPUT =================
addpath(genpath('/Users/sk3526/Library/CloudStorage/OneDrive-YaleUniversity/lab_members/Sandeep/SCaMP/main_matlab_scripts'));

folder_path = '/Volumes/Sandeep/CalciumImaging/20260528_mEXP1512_AIY_GCamP/20260528_mEXP1512_AIY_GCamP_rep_3';

plot_title = 'AFD + AIY dual-color calcium imaging';

SAVE_OUTPUT_CSV = 0;

baseline_window_s = [50 200];
baseline_percentile = 10;

ylims_afd = [];
ylims_aiy = [];

%% ================= RUN ANALYSIS =================
aligned_data = align_appended_roi_dualcolor_folder( ...
    folder_path, ...
    plot_title, ...
    SAVE_OUTPUT_CSV, ...
    baseline_window_s, ...
    baseline_percentile);

%% ================= PLOTS =================
plot_afd_aiy_2x2_summary(aligned_data, plot_title, ylims_afd, ylims_aiy);
plot_paired_afd_aiy_subplots(aligned_data, plot_title, ylims_afd, ylims_aiy);


%% ============================================================
% Main function
% ============================================================
function aligned_data = align_appended_roi_dualcolor_folder(folder_path, plot_title, SAVE_OUTPUT_CSV, baseline_window_s, baseline_percentile)

if ~isfolder(folder_path)
    error('folder_path does not exist: %s', folder_path);
end

if nargin < 2 || isempty(plot_title)
    plot_title = '';
end

if nargin < 3 || isempty(SAVE_OUTPUT_CSV)
    SAVE_OUTPUT_CSV = 0;
end

if nargin < 4 || isempty(baseline_window_s)
    baseline_window_s = [900 1100];
end

if nargin < 5 || isempty(baseline_percentile)
    baseline_percentile = 10;
end

% ------------------------------------------------------------
% Find the appended ROI file
% This can be .csv or .txt, but must contain:
% LABEL, FRAME, Mean ch1
% ------------------------------------------------------------
roi_file_path = find_appended_roi_file(folder_path);

% ------------------------------------------------------------
% Find metadata and stimulus files
% ------------------------------------------------------------
meta_files = dir(fullfile(folder_path, '*metadata.txt'));

if numel(meta_files) < 1
    error('Expected at least one *metadata.txt file in folder.');
end

meta_path = fullfile(folder_path, meta_files(1).name);
stim_path = find_stimulus_txt_file(folder_path, roi_file_path);

fprintf('\nFiles selected:\n');
fprintf('ROI data:   %s\n', roi_file_path);
fprintf('Metadata:   %s\n', meta_path);
fprintf('Stimulus:   %s\n\n', stim_path);

% ------------------------------------------------------------
% Convert appended ROI table to Frame x ROI matrix
% ------------------------------------------------------------
[roi_matrix, frame_numbers, roi_labels_by_frame] = read_appended_roi_as_matrix(roi_file_path);

nFrames = size(roi_matrix, 1);
nROI    = size(roi_matrix, 2);

fprintf('Converted ROI data to matrix: %d frames x %d ROIs\n', nFrames, nROI);

if mod(nROI, 2) ~= 0
    error('Number of ROI columns must be even. Found %d ROI columns.', nROI);
end

nPairs = nROI / 2;

% First half = AFD, second half = AIY
afd_raw = roi_matrix(:, 1:nPairs);
aiy_raw = roi_matrix(:, nPairs+1:end);

fprintf('AFD matrix: %d frames x %d ROIs\n', size(afd_raw,1), size(afd_raw,2));
fprintf('AIY matrix: %d frames x %d ROIs\n\n', size(aiy_raw,1), size(aiy_raw,2));

% ------------------------------------------------------------
% Read timing and stimulus using your existing functions
% ------------------------------------------------------------
if exist('extract_frameTimes_from_metadata', 'file') ~= 2
    error('Could not find extract_frameTimes_from_metadata on MATLAB path.');
end

if exist('extract_temps_times_from_metadata', 'file') ~= 2
    error('Could not find extract_temps_times_from_metadata on MATLAB path.');
end

if exist('align_by_index_calcium_frames_temp_1', 'file') ~= 2
    error('Could not find align_by_index_calcium_frames_temp_1 on MATLAB path.');
end

frameTimes = extract_frameTimes_from_metadata(meta_path);
Stim       = extract_temps_times_from_metadata(stim_path);

% Make sure calcium frames and metadata frame times have same length
nKeep = min([size(afd_raw,1), size(aiy_raw,1), numel(frameTimes)]);

if size(afd_raw,1) ~= numel(frameTimes)
    warning('ROI matrix has %d frames, but metadata has %d frame times. Trimming to %d frames.', ...
        size(afd_raw,1), numel(frameTimes), nKeep);
end

afd_raw = afd_raw(1:nKeep, :);
aiy_raw = aiy_raw(1:nKeep, :);
roi_matrix = roi_matrix(1:nKeep, :);
frame_numbers = frame_numbers(1:nKeep);
frameTimes = frameTimes(1:nKeep);

% ------------------------------------------------------------
% Align calcium data to temperature stimulus
% ------------------------------------------------------------
aligned_afd = align_by_index_calcium_frames_temp_1(afd_raw, frameTimes, Stim);
aligned_aiy = align_by_index_calcium_frames_temp_1(aiy_raw, frameTimes, Stim);

% ------------------------------------------------------------
% Pack output
% ------------------------------------------------------------
aligned_data = struct();

aligned_data.Time_s = aligned_afd.Time_s(:);
aligned_data.Temp   = aligned_afd.Temp(:);

aligned_data.afd_aligned = aligned_afd.afd_aligned;
aligned_data.aiy_aligned = aligned_aiy.afd_aligned;

% Kept in MATLAB only for debugging/QC.
% This is not saved as a separate CSV.
aligned_data.raw_roi_matrix = roi_matrix;
aligned_data.raw_afd_matrix = afd_raw;
aligned_data.raw_aiy_matrix = aiy_raw;
aligned_data.frame_numbers  = frame_numbers;
aligned_data.roi_labels_by_frame = roi_labels_by_frame;

aligned_data.nFrames_raw = nFrames;
aligned_data.nROI_total  = nROI;
aligned_data.nPairs      = nPairs;

aligned_data.pair_table = table( ...
    (1:nPairs)', ...
    (1:nPairs)', ...
    (nPairs+1:nROI)', ...
    'VariableNames', {'Pair', 'AFD_raw_column', 'AIY_raw_column'});

aligned_data.baseline_window_s = baseline_window_s;
aligned_data.baseline_percentile = baseline_percentile;

aligned_data.afd_dff = compute_dff_matrix( ...
    aligned_data.afd_aligned, ...
    aligned_data.Time_s, ...
    baseline_window_s, ...
    baseline_percentile);

aligned_data.aiy_dff = compute_dff_matrix( ...
    aligned_data.aiy_aligned, ...
    aligned_data.Time_s, ...
    baseline_window_s, ...
    baseline_percentile);

% ------------------------------------------------------------
% Save useful processed CSV output
% Example:
% input ROI file: 20260602_rep2.csv
% output file:    AFD_AIY_processed_data_20260602_rep2.csv
% ------------------------------------------------------------
if SAVE_OUTPUT_CSV == 1

    output_tag = get_output_tag_from_roi_file(roi_file_path);
    output_file_name = sprintf('AFD_AIY_processed_data_%s.csv', output_tag);

    aligned_out_path = fullfile(folder_path, output_file_name);
    save_aligned_output_csv(aligned_data, aligned_out_path);

    fprintf('Saved processed AFD/AIY data:\n%s\n', aligned_out_path);
end

end


%% ============================================================
% Get output tag from input ROI file name
% Example: 20260602_rep2.csv -> 20260602_rep2
% ============================================================
function output_tag = get_output_tag_from_roi_file(roi_file_path)

[~, output_tag, ~] = fileparts(roi_file_path);

% Make filename safe
output_tag = regexprep(output_tag, '[^\w]', '_');
output_tag = regexprep(output_tag, '_+', '_');
output_tag = regexprep(output_tag, '^_|_$', '');

end


%% ============================================================
% Find appended ROI file
% ============================================================
function roi_file_path = find_appended_roi_file(folder_path)

csv_files = dir(fullfile(folder_path, '*.csv'));
txt_files = dir(fullfile(folder_path, '*.txt'));

all_files = [csv_files; txt_files];

matches = {};

for k = 1:numel(all_files)

    file_name = all_files(k).name;
    file_path = fullfile(folder_path, file_name);

    if contains(file_name, 'metadata', 'IgnoreCase', true)
        continue;
    end

    if contains(file_name, 'comments', 'IgnoreCase', true)
        continue;
    end

    if contains(file_name, 'DisplaySettings', 'IgnoreCase', true)
        continue;
    end

    if contains(file_name, 'AFD_AIY_processed_data', 'IgnoreCase', true)
        continue;
    end

    if contains(file_name, 'aligned_raw_and_dff', 'IgnoreCase', true)
        continue;
    end

    if looks_like_appended_roi_file(file_path)
        matches{end+1} = file_path; %#ok<AGROW>
    end
end

if numel(matches) ~= 1
    fprintf('\nPossible ROI files found:\n');
    for i = 1:numel(matches)
        fprintf('  %s\n', matches{i});
    end
    error('Expected exactly 1 appended ROI .csv/.txt file, but found %d.', numel(matches));
end

roi_file_path = matches{1};

end


%% ============================================================
% Check file header
% ============================================================
function tf = looks_like_appended_roi_file(file_path)

tf = false;

fid = fopen(file_path, 'r');
if fid == -1
    return;
end

txt = '';
for i = 1:5
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end
    txt = [txt, ' ', line]; %#ok<AGROW>
end

fclose(fid);

tf = contains(txt, 'LABEL', 'IgnoreCase', true) && ...
     contains(txt, 'FRAME', 'IgnoreCase', true) && ...
     contains(txt, 'Mean ch1', 'IgnoreCase', true);

end


%% ============================================================
% Find stimulus txt file
% ============================================================
function stim_path = find_stimulus_txt_file(folder_path, roi_file_path)

txt_files = dir(fullfile(folder_path, '*.txt'));

stim_candidates = {};

[~, roi_name, roi_ext] = fileparts(roi_file_path);
roi_file_name = [roi_name roi_ext];

for k = 1:numel(txt_files)

    nm = txt_files(k).name;

    if strcmpi(nm, roi_file_name)
        continue;
    end

    if contains(nm, 'metadata', 'IgnoreCase', true)
        continue;
    end

    if contains(nm, 'comments', 'IgnoreCase', true)
        continue;
    end

    if contains(nm, 'DisplaySettings', 'IgnoreCase', true)
        continue;
    end

    % Prefer files such as 060226_135941_725.txt
    if ~isempty(regexp(nm, '^\d+_\d+_\d+\.txt$', 'once'))
        stim_candidates{end+1} = fullfile(folder_path, nm); %#ok<AGROW>
    end
end

if numel(stim_candidates) ~= 1
    fprintf('\nStimulus txt candidates found:\n');
    for i = 1:numel(stim_candidates)
        fprintf('  %s\n', stim_candidates{i});
    end
    error('Expected exactly 1 stimulus .txt file, but found %d.', numel(stim_candidates));
end

stim_path = stim_candidates{1};

end


%% ============================================================
% Read appended ROI table and convert to Frame x ROI matrix
% ============================================================
function [roi_matrix, frame_numbers, roi_labels_by_frame] = read_appended_roi_as_matrix(file_path)

delimiter = guess_delimiter(file_path);

opts = delimitedTextImportOptions('NumVariables', 3);
opts.DataLines = [4 Inf];
opts.Delimiter = delimiter;
opts.VariableNames = {'Label', 'Frame', 'Mean_ch1'};
opts.VariableTypes = {'string', 'double', 'double'};
opts.ExtraColumnsRule = 'ignore';
opts.EmptyLineRule = 'read';

T = readtable(file_path, opts);

% Remove empty rows
valid_rows = ~isnan(T.Frame) & ~isnan(T.Mean_ch1);
T = T(valid_rows, :);

frames = unique(T.Frame, 'stable');
nFrames = numel(frames);

counts_per_frame = zeros(nFrames, 1);

for i = 1:nFrames
    counts_per_frame(i) = sum(T.Frame == frames(i));
end

nROI = mode(counts_per_frame);

if any(counts_per_frame ~= nROI)
    warning('Not every frame has the same number of ROIs. Missing values will be filled with NaN.');
end

roi_matrix = NaN(nFrames, nROI);
roi_labels_by_frame = strings(nFrames, nROI);

for i = 1:nFrames

    idx = find(T.Frame == frames(i));
    nThis = min(numel(idx), nROI);

    roi_matrix(i, 1:nThis) = T.Mean_ch1(idx(1:nThis));
    roi_labels_by_frame(i, 1:nThis) = T.Label(idx(1:nThis));
end

frame_numbers = frames(:);

% Warn if frames are not contiguous
if any(diff(frame_numbers) ~= 1)
    warning('Frame numbers are not perfectly contiguous. Check for dropped/missing frames.');
end

end


%% ============================================================
% Guess delimiter
% ============================================================
function delimiter = guess_delimiter(file_path)

fid = fopen(file_path, 'r');
first_line = fgetl(fid);
fclose(fid);

if contains(first_line, ',')
    delimiter = ',';
elseif contains(first_line, sprintf('\t'))
    delimiter = '\t';
else
    delimiter = ',';
end

end


%% ============================================================
% Compute dF/F
% ============================================================
function dff = compute_dff_matrix(F, t, baseline_window_s, baseline_percentile)

dff = NaN(size(F));

baseline_idx = t >= baseline_window_s(1) & t <= baseline_window_s(2);

if ~any(baseline_idx)
    warning('Baseline window %.1f-%.1f s not found. Using first 10%% of recording for baseline.', ...
        baseline_window_s(1), baseline_window_s(2));

    nBase = max(1, round(0.10 * numel(t)));
    baseline_idx = false(size(t));
    baseline_idx(1:nBase) = true;
end

for c = 1:size(F, 2)

    trace = F(:, c);
    base_values = trace(baseline_idx);
    base_values = base_values(isfinite(base_values));

    if isempty(base_values)
        continue;
    end

    F0 = prctile(base_values, baseline_percentile);

    if isfinite(F0) && F0 ~= 0
        dff(:, c) = (trace - F0) ./ F0;
    end
end

end


%% ============================================================
% Save aligned output
% ============================================================
function save_aligned_output_csv(aligned_data, out_path)

t = aligned_data.Time_s(:);
Temp = aligned_data.Temp(:);

AFD_raw = aligned_data.afd_aligned;
AIY_raw = aligned_data.aiy_aligned;

AFD_dff = aligned_data.afd_dff;
AIY_dff = aligned_data.aiy_dff;

nPairs = min(size(AFD_raw, 2), size(AIY_raw, 2));

AFD_raw = AFD_raw(:, 1:nPairs);
AIY_raw = AIY_raw(:, 1:nPairs);
AFD_dff = AFD_dff(:, 1:nPairs);
AIY_dff = AIY_dff(:, 1:nPairs);

outMat = [t, Temp, AFD_raw, AIY_raw, AFD_dff, AIY_dff];

varNames = cell(1, size(outMat, 2));
varNames{1} = 'Time_s';
varNames{2} = 'Temp';

col = 3;

for i = 1:nPairs
    varNames{col} = sprintf('AFD_%02d_raw', i);
    col = col + 1;
end

for i = 1:nPairs
    varNames{col} = sprintf('AIY_%02d_raw', i);
    col = col + 1;
end

for i = 1:nPairs
    varNames{col} = sprintf('AFD_%02d_dFF', i);
    col = col + 1;
end

for i = 1:nPairs
    varNames{col} = sprintf('AIY_%02d_dFF', i);
    col = col + 1;
end

outT = array2table(outMat, 'VariableNames', varNames);
writetable(outT, out_path);

end


%% ============================================================
% Plot AFD and AIY summary as one 2x2 figure
% Top row: AFD traces + AFD heatmap
% Bottom row: AIY traces + AIY heatmap
% ============================================================
function plot_afd_aiy_2x2_summary(aligned_data, plot_title, ylims_afd, ylims_aiy)

t = aligned_data.Time_s(:);
Temp = aligned_data.Temp(:);

AFD = aligned_data.afd_dff;
AIY = aligned_data.aiy_dff;

xMax = max(t(isfinite(t)));

AFD_TRUE_COLOR  = [1 0 0];
AFD_TRACE_COLOR = [1.0 0.4 0.4];
AFD_MEAN_COLOR  = [0.45 0 0];

AIY_TRUE_COLOR  = [0 0.60 0];
AIY_TRACE_COLOR = [0.4 0.76 0.4];
AIY_MEAN_COLOR  = [0 0.30 0];

TEMP_COLOR = [0 0 0];

figure('Color', 'w', 'Name', 'AFD_AIY_2x2_summary');
tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

sgtitle(tl, plot_title, 'FontWeight', 'normal', 'Interpreter', 'none');

% ---------------- AFD traces ----------------
ax1 = nexttile(tl, 1);
plot_single_channel_traces_with_temp( ...
    ax1, t, Temp, AFD, ...
    AFD_TRUE_COLOR, AFD_TRACE_COLOR, AFD_MEAN_COLOR, TEMP_COLOR, ...
    'AFD together', ...
    'AFD \DeltaF/F_0', ...
    'Individual AFDs', ...
    ylims_afd, ...
    xMax);

% ---------------- AFD heatmap ----------------
ax2 = nexttile(tl, 2);
plot_single_channel_heatmap( ...
    ax2, t, AFD, ...
    'AFD heatmap', ...
    'AFD', ...
    xMax);

% ---------------- AIY traces ----------------
ax3 = nexttile(tl, 3);
plot_single_channel_traces_with_temp( ...
    ax3, t, Temp, AIY, ...
    AIY_TRUE_COLOR, AIY_TRACE_COLOR, AIY_MEAN_COLOR, TEMP_COLOR, ...
    'AIY together', ...
    'AIY \DeltaF/F_0', ...
    'Individual AIYs', ...
    ylims_aiy, ...
    xMax);

% ---------------- AIY heatmap ----------------
ax4 = nexttile(tl, 4);
plot_single_channel_heatmap( ...
    ax4, t, AIY, ...
    'AIY heatmap', ...
    'AIY', ...
    xMax);

end


%% ============================================================
% Helper: traces + mean + temperature
% ============================================================
function plot_single_channel_traces_with_temp(ax, t, Temp, X, true_trace_color, trace_color, mean_color, TEMP_COLOR, panel_title, y_label, individual_label, ylimits, xMax)

hold(ax, 'on');

mean_trace = mean(X, 2, 'omitnan');

yyaxis(ax, 'left');

% Individual traces
for i = 1:size(X, 2)
    plot(ax, t, X(:, i), '-', 'Color', trace_color, 'LineWidth', 1.0);
end

% Legend proxy for individual traces
p_ind = plot(ax, nan, nan, '-', 'Color', true_trace_color, 'LineWidth', 1.0);

% Mean trace
p_mean = plot(ax, t, mean_trace, '-', 'Color', mean_color, 'LineWidth', 2);

ylabel(ax, y_label, 'Color', true_trace_color);

if ~isempty(ylimits)
    ylim(ax, ylimits);
end

yyaxis(ax, 'right');
p_temp = plot(ax, t, Temp, 'k-', 'LineWidth', 1.5);
ylabel(ax, 'Temperature', 'Color', TEMP_COLOR);

xlabel(ax, 'Time (s)');
title(ax, panel_title, 'FontWeight', 'normal', 'Interpreter', 'none');

xlim(ax, [0 xMax]);

ax.YAxis(1).Color = true_trace_color;
ax.YAxis(2).Color = TEMP_COLOR;

box(ax, 'off');
set(ax, 'FontSize', 14);

legend(ax, [p_ind p_mean p_temp], {individual_label, 'Mean', 'Temperature'}, 'Location', 'northwest');
legend(ax, 'boxoff');

end


%% ============================================================
% Helper: heatmap
% ============================================================
function plot_single_channel_heatmap(ax, t, X, panel_title, y_label, xMax)

imagesc(ax, t, 1:size(X, 2), X');
set(ax, 'YDir', 'normal');

xlabel(ax, 'Time (s)');
ylabel(ax, y_label);
title(ax, panel_title, 'FontWeight', 'normal', 'Interpreter', 'none');

xlim(ax, [0 xMax]);

colormap(ax, get_viridis_colormap(256));

cb = colorbar(ax);
cb.Label.String = '\DeltaF/F_0';

box(ax, 'off');
set(ax, 'FontSize', 14);

end


%% ============================================================
% Viridis-like colormap
% This avoids errors if MATLAB does not have viridis built in
% ============================================================
function cmap = get_viridis_colormap(n)

if nargin < 1
    n = 256;
end

base = [
    0.2670 0.0049 0.3294
    0.2823 0.1409 0.4575
    0.2539 0.2653 0.5299
    0.2068 0.3718 0.5531
    0.1636 0.4711 0.5581
    0.1276 0.5669 0.5506
    0.1347 0.6586 0.5176
    0.2669 0.7488 0.4406
    0.4775 0.8214 0.3182
    0.7414 0.8734 0.1496
    0.9932 0.9062 0.1439
];

x = linspace(0, 1, size(base, 1));
xi = linspace(0, 1, n);

cmap = interp1(x, base, xi, 'linear');

end


%% ============================================================
% Plot paired AFD and AIY in subplots
% 2 rows, optimized number of columns
% ============================================================
function plot_paired_afd_aiy_subplots(aligned_data, plot_title, ylims_afd, ylims_aiy)

t = aligned_data.Time_s(:);

AFD = aligned_data.afd_dff;
AIY = aligned_data.aiy_dff;

nPairs = min(size(AFD, 2), size(AIY, 2));

AFD_COLOR = [1 0 0];
AIY_COLOR = [0 0.60 0];

xMax = max(t(isfinite(t)));
xlims = [0 xMax];

nRows = 2;
nCols = ceil(nPairs / nRows);

figure('Color', 'w', 'Name', 'AFD_AIY_paired_subplots');
tl = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

sgtitle(tl, 'Paired AFD and AIY', ...
    'FontWeight', 'normal', ...
    'Interpreter', 'none');

for p = 1:nPairs

    ax = nexttile(tl, p);
    hold(ax, 'on');

    yyaxis(ax, 'left');
    pA = plot(ax, t, AFD(:, p), '-', 'Color', AFD_COLOR, 'LineWidth', 1.5);
    ylabel(ax, 'AFD \DeltaF/F_0', 'Color', AFD_COLOR);

    if nargin >= 3 && ~isempty(ylims_afd)
        ylim(ax, ylims_afd);
    end

    yyaxis(ax, 'right');
    pI = plot(ax, t, AIY(:, p), '-', 'Color', AIY_COLOR, 'LineWidth', 1.5);
    ylabel(ax, 'AIY \DeltaF/F_0', 'Color', AIY_COLOR);

    if nargin >= 4 && ~isempty(ylims_aiy)
        ylim(ax, ylims_aiy);
    end

    xlim(ax, xlims);
    xlabel(ax, 'Time (s)');

    ax.YAxis(1).Color = AFD_COLOR;
    ax.YAxis(2).Color = AIY_COLOR;

    box(ax, 'off');
    set(ax, 'FontSize', 12);

    if p == 1
        legend(ax, [pA pI], {'AFD', 'AIY'}, 'Location', 'northwest');
        legend(ax, 'boxoff');
    end
end

% Hide unused tiles
for k = nPairs+1:nRows*nCols
    ax = nexttile(tl, k);
    axis(ax, 'off');
end

end
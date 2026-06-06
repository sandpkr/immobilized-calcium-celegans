function aligned_data = align_by_index_calcium_frames_temp_1(afd_data, frameTimes, Stim)
% align_by_index (3 inputs)
% Aligns imaging calcium traces to stimulus time.
%
% Inputs:
%   afd_data   : [nFrames x nWorms] (or [nFrames x nCols] where each col is a worm/ROI)
%   frameTimes : [nFrames x 1] absolute time in ms for each imaging frame
%   Stim       : struct with fields:
%                  Stim.Time [nStim x 1] absolute time in ms
%                  Stim.Temp [nStim x 1] temperature
%
% Output aligned_data:
%   aligned_data.afd_aligned : [nStim x nWorms] calcium aligned to Stim.Time
%   aligned_data.Time_ms     : [nStim x 1] stimulus absolute time (ms)
%   aligned_data.Time_s      : [nStim x 1] time relative to start (s)
%   aligned_data.Temp        : [nStim x 1] temperature at each stimulus sample
%   aligned_data.frame_idx   : [nStim x 1] matched imaging frame index for each stim sample

% Ensure column vectors
frameTimes = frameTimes(:);
stimTime   = Stim.Time(:);
stimTemp   = Stim.Temp(:);

if size(afd_data,1) ~= numel(frameTimes)
    error('afd_data rows (%d) must match frameTimes length (%d).', size(afd_data,1), numel(frameTimes));
end
if numel(stimTime) ~= numel(stimTemp)
    error('Stim.Time and Stim.Temp must be same length.');
end

nFrames = numel(frameTimes);
nStim   = numel(stimTime);
nWorms  = size(afd_data,2);

afd_aligned = nan(nStim, nWorms);
frame_idx   = nan(nStim, 1);

% For each stim time, pick LAST imaging frame at or before stim time.
for i = 1:nStim
    t0 = stimTime(i);
    j = find(frameTimes <= t0, 1, 'last');

    if isempty(j)
        j = 1;
    end
    if j > nFrames
        j = nFrames;
    end

    afd_aligned(i,:) = afd_data(j,:);
    frame_idx(i) = j;
end

aligned_data = struct();
aligned_data.afd_aligned = afd_aligned;
aligned_data.Time_ms     = stimTime;
aligned_data.Time_s      = (stimTime - stimTime(1)) / 1000;
aligned_data.Temp        = stimTemp;
aligned_data.frame_idx   = frame_idx;

end

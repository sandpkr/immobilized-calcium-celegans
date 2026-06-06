function Stim_data = extract_temps_times_from_metadata(stim_filename)
    % Initializes a struct to hold the temperature and time data.
    Stim_data = struct('Temp', [], 'Time', []);
    
    % Opens the stimulus file for reading.
    fid = fopen(stim_filename, 'r');
    
    % Extracts just the filename from the full path.
    path_parts = strsplit(stim_filename, '/');
    filename_only = path_parts{end};
    
    % Parses the filename to get the starting hour and minute.
    rep_to_time = strsplit(strrep(filename_only, '.txt', ''), '_');
    hour = str2double(rep_to_time{2}(1:2));
    minute = str2double(rep_to_time{2}(3:4));
    
    % Converts the starting time to milliseconds.
    hour_min = (hour * 3600 * 1000) + (minute * 60 * 1000);
    
    % Initializes variables to track elapsed time and handle minute rollovers.
    elapsed_time_check = 0;
    elapsed_min = 0;
    
    % Reads the file line by line.
    tline = fgetl(fid);
    while ischar(tline)
        % Splits the line by tabs.
        pieces = strsplit(tline, '\t');
        
        % Skips lines that don't have the expected format.
        if length(pieces) < 2 || length(pieces) > 4
            tline = fgetl(fid);
            continue;
        end
        
        % Appends the temperature value.
        Stim_data.Temp(end+1) = str2double(pieces{1});
        
        % Parses the seconds and milliseconds from the timestamp.
        time_parts = strsplit(pieces{4}, '.');
        sec = str2double(time_parts{1});
        ms = str2double(time_parts{2});
        
        % Calculates the elapsed time in milliseconds for the current line.
        elapsed_time = (elapsed_min * 60 * 1000) + (sec * 1000) + ms;
        
        % This logic handles the minute rollover. If the current elapsed time
        % is less than the previous one, it means a minute has passed.
        if elapsed_time > elapsed_time_check
            Stim_data.Time(end+1) = hour_min + elapsed_time;
            elapsed_time_check = elapsed_time;
        else
            elapsed_min = elapsed_min + 1;
            elapsed_time = (elapsed_min * 60 * 1000) + (sec * 1000) + ms;
            Stim_data.Time(end+1) = hour_min + elapsed_time;
            elapsed_time_check = elapsed_time;
        end
        
        % Reads the next line.
        tline = fgetl(fid);
    end
    
    % Closes the file.
    fclose(fid);
end
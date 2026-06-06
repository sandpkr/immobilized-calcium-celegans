function frameTimes = extract_frameTimes_from_metadata(metadata_filename)
    frameTimes = [];

    fid = fopen(metadata_filename);
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    metadata = jsondecode(str);

    if str2double(metadata.Summary.MicroManagerVersion(1)) == 1
        % % disp('MicroManager v1');
        hour = str2double(metadata.FrameKey_0_0_0.Time(12:13));
        minute = str2double(metadata.FrameKey_0_0_0.Time(15:16));
        second = str2double(metadata.FrameKey_0_0_0.Time(18:19));
        start_time = (hour*60*60*1000) + (minute*60*1000) + (second*1000) - metadata.FrameKey_0_0_0.ElapsedTime_ms;
    elseif str2double(metadata.Summary.MicroManagerVersion(1)) > 1
        % % disp('MicroManager v2 or greater');
        hour = str2double(metadata.Summary.StartTime(12:13));
        minute = str2double(metadata.Summary.StartTime(15:16));
        second = str2double(metadata.Summary.StartTime(18:19));
        start_time = (hour*60*60*1000) + (minute*60*1000) + (second*1000) + str2double(metadata.Summary.StartTime(21:23));
    end

    fields = fieldnames(metadata);
    for i = 1:numel(fields)
        if startsWith(fields{i}, "FrameKey")
            frame_name = fields{i};
            frameTimes = [frameTimes; start_time + metadata.(frame_name).ElapsedTime_ms];
        end
    end
end
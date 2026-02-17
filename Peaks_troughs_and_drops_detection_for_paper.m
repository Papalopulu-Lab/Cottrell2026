% Peak and Trough detection for amplitude and fold change in Fig 4

Array_of_cells = % Here is the data

Big_Trend = cell(1,numel(Array_of_cells));
Detrended_array_of_cells = cell(1,numel(Array_of_cells));
Ultradian_trend = cell(1,numel(Array_of_cells));
Peaks_of_trend = cell(1,numel(Array_of_cells));
Peak_locations_of_trace = cell(1,numel(Array_of_cells));
Peak_widths_of_trace = cell(1,numel(Array_of_cells));
Peaks_of_trace = cell(1,numel(Array_of_cells));
Troughs_of_trend = cell(1,numel(Array_of_cells));
Trough_locations_of_trace = cell(1,numel(Array_of_cells));
Trough_widths_of_trace = cell(1,numel(Array_of_cells));
Troughs_of_trace = cell(1,numel(Array_of_cells));
Peaks_and_troughs = cell(1,numel(Array_of_cells));
Peak_to_trough_difference = cell(1,numel(Array_of_cells));
Peaks_and_troughs_locs = cell(1,numel(Array_of_cells));

Smoothing_window_Big = % Larger value as mentioned in paper
Smoothing_window_Ultradian= % Smaller value as mentioned in paper

for Cell_index = 1:numel(Array_of_cells) % loop through every trace

    % Find the circadian-like trend

    Big_Trend{Cell_index} = smoothdata(Array_of_cells{Cell_index},'gaussian',Smoothing_window_Big);

    % We can now find peaks at the circadian level with the smoothened trace

    % Detrend the data

    Detrended_array_of_cells{Cell_index} = Array_of_cells{Cell_index}-Big_Trend{Cell_index};


    % We use findpeaks with a minimum prominance threshold to determine
    % where our circadian peaks lie

    [Peaks_of_trend{Cell_index},Peak_locations_of_trace{Cell_index},Peak_widths_of_trace{Cell_index}] = ...
        findpeaks(Big_Trend{Cell_index},'MinPeakProminence',0.5);
    [Troughs_of_trend{Cell_index},Trough_locations_of_trace{Cell_index},Trough_widths_of_trace{Cell_index}] = ...
        findpeaks(-Big_Trend{Cell_index},'MinPeakProminence',0.5);


    % We assign the peaks to be the raw values at the time of the peak
    % found in the artifiical smoothened trace

    Peaks_of_trace{Cell_index} = Array_of_cells{Cell_index}(Peak_locations_of_trace{Cell_index});
    Troughs_of_trace{Cell_index} = Array_of_cells{Cell_index}(Trough_locations_of_trace{Cell_index});


    % Below is code to ensure our peaks and troughs are ordered in the
    % correct way: i.e. we determine peak to trough difference per cell by
    % comparing the height at a peak subtracted by the subsequent trough.
    % This means that we must make adjustments if our code does not find a
    % peak before a trough: if this happens, we consider the firs tpeak of
    % the trace to be the highest value of the trace before the first
    % torugh.


    if ~isempty(Peaks_of_trace{Cell_index}) == 1 & ~isempty(Troughs_of_trace{Cell_index}) == 1

        if Peak_locations_of_trace{Cell_index}(1) < Trough_locations_of_trace{Cell_index}(1)

            Peaks_and_troughs{Cell_index}(1:length(Peaks_of_trace{Cell_index}),1) =...
                Peaks_of_trace{Cell_index};

            Peaks_and_troughs{Cell_index}(1:length(Troughs_of_trace{Cell_index}),2) =...
                Troughs_of_trace{Cell_index};

            Peaks_and_troughs_locs{Cell_index}(1:length(Peak_locations_of_trace{Cell_index}),1) =...
                Peak_locations_of_trace{Cell_index};

            Peaks_and_troughs_locs{Cell_index}(1:length(Trough_locations_of_trace{Cell_index}),2) =...
                Trough_locations_of_trace{Cell_index};

        else

            Peaks_and_troughs{Cell_index}(1,1) =...
                max(Array_of_cells{Cell_index}(1:Trough_locations_of_trace{Cell_index}(1)));

            Peaks_and_troughs{Cell_index}(2:1+length(Peaks_of_trace{Cell_index}),1) =...
                Peaks_of_trace{Cell_index};

            Peaks_and_troughs{Cell_index}(1:length(Troughs_of_trace{Cell_index}),2) =...
                Troughs_of_trace{Cell_index};

            Peaks_and_troughs_locs{Cell_index}(1,1) =...
                find(Array_of_cells{Cell_index}==...
                max(Array_of_cells{Cell_index}(1:Trough_locations_of_trace{Cell_index}(1))),1);

            Peaks_and_troughs_locs{Cell_index}(2:1+length(Peaks_of_trace{Cell_index}),1) =...
                Peak_locations_of_trace{Cell_index};

            Peaks_and_troughs_locs{Cell_index}(1:length(Trough_locations_of_trace{Cell_index}),2) =...
                Trough_locations_of_trace{Cell_index};

        end


    elseif isempty(Peaks_of_trace{Cell_index}) == 1 & isscalar(Troughs_of_trace{Cell_index})

        Peaks_and_troughs{Cell_index}(1,1) =...
            max(Array_of_cells{Cell_index}(1:Trough_locations_of_trace{Cell_index}(1)));

        Peaks_and_troughs{Cell_index}(1:length(Troughs_of_trace{Cell_index}),2) =...
            Troughs_of_trace{Cell_index};

        Peaks_and_troughs_locs{Cell_index}(1,1) =...
            find(Array_of_cells{Cell_index}==...
            max(Array_of_cells{Cell_index}(1:Trough_locations_of_trace{Cell_index}(1))));

        Peaks_and_troughs_locs{Cell_index}(1:length(Trough_locations_of_trace{Cell_index}),2) =...
            Trough_locations_of_trace{Cell_index};

    end

    if ~isempty(Peaks_of_trace{Cell_index}) == 1

        Peaks_and_troughs{Cell_index}(Peaks_and_troughs{Cell_index} == 0) = ...
            min(Array_of_cells{Cell_index}(Peak_locations_of_trace{Cell_index}(end,1):end));

        Peaks_and_troughs_locs{Cell_index}(Peaks_and_troughs_locs{Cell_index} == 0) = ...
            find(Array_of_cells{Cell_index}==...
            min(Array_of_cells{Cell_index}(Peak_locations_of_trace{Cell_index}(end,1):end)),1,'last');

    end

    if ~isempty(Peaks_and_troughs{Cell_index}) == 1

        % Absolute value difference (amplitude)

        Peak_to_trough_difference{Cell_index} = Peaks_and_troughs{Cell_index}(:,1)-Peaks_and_troughs{Cell_index}(:,2);
    end



end

% We smooth again to be able to find ultradian peaks and repeat the above process for the ultradian trend

for Cell_index = 1:numel(Array_of_cells) % loop through every trace

    Ultradian_trend{Cell_index} = smoothdata(smoothdata(Detrended_array_of_cells{Cell_index},'gaussian',Smoothing_window_Ultradian),'gaussian',Smoothing_window_Ultradian);

    % Here would go the repeated code as above but inputting
    % Ultradian_trend into findpeaks

end


%% CDK2 and Hes1 peak and trough detection for Fig 3

Array_of_cells = % Hes1 or CDK2 data goes here
Time_of_max_peak = nan(1,size(Data,2));
Time_of_min_trough = nan(1,size(Data,2));

for Cell_index = 1:numel(Array_of_cells) % loop through every trace

    [pks, locs] = findpeaks(Array_of_cells{Cell_index}(15:end-15));
    if ~isempty(pks)
        Time_of_max_peak(Cell_index) = locs(pks==max(pks));
    end
    [troughs, trough_locs] = findpeaks(-Array_of_cells{Cell_index}(15:end-15));

    if ~isempty(troughs)
        Time_of_min_trough(Cell_index) = trough_locs(troughs==max(troughs));
    end

end

%% Hes1 dip and P21 drop code for Fig 3

% P21 drop identification

Hes1_raw_traces = % Here we put our Hes1 data;
P21_traces_array = % Here we put our Hes1 data;
% these are from the same cell so can go through the same loop


% Preallocation
Hes1_smoothened_traces = cell(1,length(Hes1_raw_traces));
Turning_points = cell(1,length(Hes1_raw_traces));
Detected_dip = zeros(1,length(Hes1_raw_traces));
Time_from_start_of_cell_cycle_until_dip = zeros(1,length(Hes1_raw_traces));
Time_from_dip_until_end_of_cell_cyle = zeros(1,length(Hes1_raw_traces));

P21_scaled_traces_array_zscore = cell(1,length(Hes1_raw_traces));

P21_drop_smooth = cell(1,length(Hes1_raw_traces));

P21_drop_time = cell(1,length(Hes1_raw_traces));

% Data observed every 20 mins
time_conversion_scalar = 3;

% We require the correct smoothing to successfully capture the drop
P21_drop_smoothing_window=7;

% Loop through all cells
for cell_index = 1:numel(Hes1_raw_traces)

% Normalise cells
    P21_scaled_traces_array_zscore{cell_index} = zscore(P21_traces_array{cell_index});

% Smooth data 
    P21_drop_smooth{cell_index} = smoothdata(P21_scaled_traces_array_zscore{cell_index},...
        'gaussian',P21_drop_smoothing_window);

% Determine where the drop happens. Defined as the point at which the
% difference between subsequent points in the smoothened trace is
% mninimised. However we only chhose this point if its within the first
% 3/4's of the trace length to avoid artefacts. 

    P21_drop_time(cell_index) = ...
        find(diff(P21_drop_smooth{cell_index}(1:round(3*numel(P21_drop_smooth{cell_index})/4)))...
        ==min(diff(P21_drop_smooth{cell_index}(1:round(3*numel(P21_drop_smooth{cell_index})/4)))))/time_conversion_scalar;


    % Hes1 dip identification


    % Set a moving frame length for the Savitzgy-Golay filter, we choose
    % this to be at a quarter of the length of each trace. The fram length
    % must be odd in order to work, hence the if loop. We also set the
    % order of polynmomial as 3, this means a local cubic fit is performed
    % for each time point as the frame moves across the trace.

    if mod(round((1/4)*(length(Hes1_raw_traces{Cell_index}))),2) == 1
        % if a quarter of the trace length is odd, then use this value
        frame_length_for_filter = round((1/4)*(length(Hes1_raw_traces{Cell_index})));
    else
        % if a quarter of the trace length is even, then use one less than this value
        frame_length_for_filter = round((1/4)*(length(Hes1_raw_traces{Cell_index})))-1; %MUST BE ODD!
    end

    % Apply the Savitzgy-Golay filter

    Hes1_smoothened_traces{Cell_index} = sgolayfilt(Hes1_raw_traces{Cell_index},3,frame_length_for_filter);

    % Chop off the ends of the smoothened trace, to abandon smoothing artefacts

    Truncated_smoothened_trace = Hes1_smoothened_traces{Cell_index}(3:end-2);

    Turning_points{Cell_index} = [];

    % Find the turning points of the smooth trace. This is done by
    % multiplying the difference of three consecutive time points.
    % For instance, if we have an upward trend and a downward trend either
    % side of a time point then this product will be negative, i.e. a
    % turning point. We record every time this is the case.

    for time_index = 2:length(Truncated_smoothened_trace)-2
        if (Truncated_smoothened_trace(time_index) - Truncated_smoothened_trace(time_index+1))*...
                (Truncated_smoothened_trace(time_index+1) - Truncated_smoothened_trace(time_index+2)) < 0
            Turning_points{Cell_index} = [Turning_points{Cell_index},(time_index+3)];
        end
    end

    turning_points_not_at_start_of_cell_cycle = [];

    % Some turning points are captured early in the cell cycle, we know that
    % cells sometimes have an early dip in expression, yet this not the
    % only dip in the whole trace and not the one we are looking for. To
    % account for this, we disgard turning points in the first 30% of the
    % cell cycle.

    for t_point = 1:length(Turning_points{Cell_index})
        if  Turning_points{Cell_index}(t_point) > (length(Hes1_smoothened_traces{Cell_index}))*(2/10) &&...
                Turning_points{Cell_index}(t_point) < (length(Hes1_smoothened_traces{Cell_index}))*(7/10)
            turning_points_not_at_start_of_cell_cycle = [turning_points_not_at_start_of_cell_cycle, Turning_points{Cell_index}(t_point)];
        end
    end

    % Our detected dip is now considered to be the lowest (minimum) of such remaining
    % turning points.

    Detected_dip(Cell_index) = find(Hes1_smoothened_traces{Cell_index} == min(Hes1_smoothened_traces{Cell_index}(turning_points_not_at_start_of_cell_cycle))); %
    Time_from_start_of_cell_cycle_until_dip(Cell_index) = Detected_dip(Cell_index)/time_conversion_scalar;

end


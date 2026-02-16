% Within this section, we use Lomb-Scargle periodogram to estimate the
% dominant periodicity of each time series trace

time_conversion_scalar = 3; % Time points of our data was every 20 mins
max_length = 300; % Arbitrary large value for interpolation purposes

% Pre-Allocate

LSP_power = cell(1,size(Cells,2));
LSP_frequency = cell(1,size(Cells,2));
Dominant_period = zeros(1,size(Cells,2));
Dominant_power = zeros(1,size(Cells,2));
Dominant_circadian_power = zeros(1,size(Cells,2));
VoS_freq = cell(1,size(Cells,2));
VOS_matrix = zeros(max_length,size(Cells,2));
LSP_power_interpolated  = zeros(max_length,size(Cells,2));
Normalised_data= nan(max_length,size(Cells,2));


% Loop through cells

for Cell_index = 1:size(Cells,2)


    y = % Hear we put our data

    time_vector = (1:length(y))/time_conversion_scalar;

    freq_vec = 0:0.001:0.5;

    [LSP_power{Cell_index}, LSP_frequency{Cell_index}] = plomb(y,time_vector,freq_vec,'normalized');

    Dominant_frequency(Cell_index) = LSP_frequency{Cell_index}(LSP_power{Cell_index}==max(LSP_power{Cell_index}));

    Dominant_period(Cell_index) = (1/Dominant_frequency(Cell_index))/time_conversion_scalar;

    Dominant_power(Cell_index) =...
        max(LSP_power{Cell_index});

    % interpolate power and freq in order to average correctly

    VoS_freq{Cell_index}=linspace(LSP_frequency{Cell_index}(1),(LSP_frequency{Cell_index}(end),max_length);

    VOS_matrix(:,Cell_index) =VoS_freq{Cell_index};

    LSP_power_interpolated(:,Cell_index) =interp1(LSP_frequency{Cell_index},LSP_power{Cell_index},VoS_freq{Cell_index});


end


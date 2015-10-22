clear;

%{
Calculates sum of squared pressures over time for all sensor grid points.
This is proportional to the intensity field.
%}

fprintf('Starting intensity calculation');

h5input = 'kwave_input_data.h5';
h5output = 'kwave_output_data.h5';

%% Set transducer excitation frequency here.
% * Need a better way to do this so we don't have to hard-code it each time
f0 = 2.36e6; % Hz

%% Get frequency coordinates from sampling freq (fs) limits.
Nt = h5read(h5input, '/Nt');
dt = h5read(h5input, '/dt'); % s
fs = 1/dt; % Hz
f = linspace(0, fs, Nt); % Hz

% Only need up to the 2nd harmonic. Add some overhead to make sure we fully
% capture all 3 frequencies.
max_harmonic = 2;
max_freq = (max_harmonic+1.5)*f0; % Hz
max_f_index = find(f<max_freq, 1, 'last')

sensor_mask_indices = h5read(h5input, '/sensor_mask_index');
num_sensor_points = length(sensor_mask_indices);

p = zeros(Nt, num_sensor_points);

%% Read in entire pressure data at each sensor
for timestep=double(1:Nt)
    p(timestep, :) = h5read(h5output, '/p', [1, timestep, 1], ...
                            [num_sensor_points, 1, 1]);
    % Progress bar.
    timestep
    %if mod(sensor_index, round(length(intensity)/100)) == 0
    %    fprintf('%.2f%% complete.\n', 100*timestep/length(intensity))
    %end
end

%% Calculate intensity field and FFT (wrt time) at each sensor location.
% Note: need to write something in Python and see how much faster it is.
intensity = zeros(num_sensor_points, 1);
p_fft = zeros(num_sensor_points, max_f_index);
i = 1;

for p_sensor=p    
    intensity(i) = sumsqr(p_sensor);
    
    % Do FFT, make sure length of the pressure vector is a power of 2 to
    % optimize the FFT.
    p_sensor_fft = fft(double(p_sensor), 2^nextpow2(double(Nt)));
    p_fft(i, :) = p_sensor_fft(1:max_f_index);
    
    i = i + 1;
end

%% Map intensities to x-z and x-y planes based on sensor mask indices
Nx = h5read(h5input,'/Nx');
Ny = h5read(h5input,'/Ny');
Nz = h5read(h5input,'/Nz');

% lateral (x-y) plane
lat_plane_intensities = zeros(Nx, Ny); 
lat_plane_fft = zeros(Nx, Ny, max_f_index);
% elevational (x-z) plane
ele_plane_intensities = zeros(Nx, Nz);
ele_plane_fft = zeros(Nx, Nz, max_f_index);

[ax, lat, ele] = ind2sub([Nx, Ny, Nz], sensor_mask_indices);

for i=1:length(intensity)
    if mod(i, round(length(intensity)/100)) == 0
        fprintf('%.2f%% complete\n', 100*i/length(ax));
    end
    x = ax(i);
    y = lat(i);
    z = ele(i);
    
    % center elevational position means we're on lateral plane
    if z == Nz/2
        lat_plane_intensities(x, y) = intensity(i);
        lat_plane_fft(x, y, :) = p_fft(i, :);
    end

    % center lateral position means we're on elevational plane
    if y == Ny/2
        ele_plane_intensities(x, z) = intensity(i);
        ele_plane_fft(x, z, :) = p_fft(i, :);
    end
    
    if z ~= Nz/2 && y ~= Ny/2
        sprintf('fail');
        break
    end
end


save('kwave_fft_intensity_data.mat', 'Nx', 'Ny', 'Nz', 'ele_plane_fft', 'ele_plane_intensities', 'lat_plane_fft', 'lat_plane_intensities', 'p_fft', 'intensity')

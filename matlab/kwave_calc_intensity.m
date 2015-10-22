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

sensor_mask_indices = h5read(h5input, '/sensor_mask_index');
num_sensor_points = length(sensor_mask_indices);

intensity = zeros(num_sensor_points, 1);

for timestep=double(1:Nt)
    p_timestep = h5read(h5output, '/p', [1, timestep, 1], ...
                            [num_sensor_points, 1, 1]);
    intensity = intensity + p_timestep.^2;
    timestep
    % Progress bar.
    %if mod(sensor_index, round(length(intensity)/100)) == 0
    %    fprintf('%.2f%% complete.\n', 100*timestep/length(intensity))
    %end
end
    
save('kwave_intensity_data_raw.mat', 'intensity')

%% Map intensities to x-z and x-y planes based on sensor mask indices
Nx = h5read(h5input,'/Nx');
Ny = h5read(h5input,'/Ny');
Nz = h5read(h5input,'/Nz');

% lateral (x-y) plane
lat_plane_intensities = zeros(Nx, Ny); 
% elevational (x-z) plane
ele_plane_intensities = zeros(Nx, Nz);

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
    end

    % center lateral position means we're on elevational plane
    if y == Ny/2
        ele_plane_intensities(x, z) = intensity(i);
    end
    
    if z ~= Nz/2 && y ~= Ny/2
        sprintf('fail');
        break
    end
end

save('kwave_intensity_data_full.mat', 'Nx', 'Ny', 'Nz', 'ele_plane_intensities', 'lat_plane_intensities')

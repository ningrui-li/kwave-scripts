% ch41_80pwr_txBeam, based on Simulating Ultrasound Beam Patterns Example:
% example_us_beam_patterns
%
% Original comments are maintained below, and I have edited along the way to
% simulate the CH4-1 using source pressures that Ned Rouze measured in:
% /luscinia/ncr2/intensity_measurements/20121130_ch4-1_80pwr_focus_face_scans_sonic/
% face_scan_80pwr_sonic_hydrophone/p0p0mm/ 
% 
% Mark Palmeri
% 2015-05-11
%
% This example shows how the nonlinear beam pattern from an ultrasound
% transducer can be modelled. It builds on the Defining An Ultrasound
% Transducer and Simulating Transducer Field Patterns examples. 
%
% author: Bradley Treeby
% date: 27th July 2011
% last update: 25th September 2012
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>. 

clear all;
% k-Wave should be added to path in startup.m
% addpath('/home/mlp6/matlab/k-Wave')
addpath /luscinia/nl91/kwave/k-Wave/matlab/k-Wave/

h5input = 'kwave_input_data.h5';
h5output = 'kwave_output_data.h5';

% simulation settings
DATA_CAST = 'single';           % set to 'single' or 'gpuArray-single' to speed up computations
MASK_PLANE = 'xy-xz';           % set to 'xy' or 'xz' to generate the beam pattern in different planes
USE_STATISTICS = false;         % set to true to compute the rms or peak beam patterns, set to false 
                                % to compute the harmonic beam patterns

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
PML_X_SIZE = 20;            % [grid points]
PML_Y_SIZE = 10;            % [grid points]
PML_Z_SIZE = 10;            % [grid points]

% set total number of grid points (not including the PML)
Nx = 2048 - 2*PML_X_SIZE;    % [grid points]
Ny = 2048 - 2*PML_Y_SIZE;    % [grid points]
Nz = 512 - 2*PML_Z_SIZE;     % [grid points]

% set desired grid size in the x-direction (not including the PML)
x = 59e-3;                  % [m]

% calculate the spacing between the grid points
dx = x/Nx;                  % [m]
dy = dx;                    % [m]
dz = dx;                    % [m]

% create the k-space grid
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% define the properties of the propagation medium
medium.sound_speed = 1540;      % [m/s]
medium.density = 1030;          % [kg/m^3]
medium.alpha_coeff = 0.30;      % [dB/(MHz^y cm)]
medium.alpha_power = 1.1;
medium.BonA = 5;

% create the time array
kgrid.t_array = makeTime(kgrid, medium.sound_speed);
%t_end = 50e-6;                  % [s]
%kgrid.t_array = makeTime(kgrid, medium.sound_speed, [], t_end);

% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================

% define properties of the input signal
source_strength = 0.4e6;        % [Pa]
tone_burst_freq = 2.5e6;    	% [Hz]
tone_burst_cycles = 7;

% create the input signal using toneBurst 
input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles, 'Envelope', 'Rectangular');

% scale the source magnitude by the source_strength divided by the
% impedance (the source is assigned to the particle velocity)
input_signal = (source_strength./(medium.sound_speed*medium.density)).*input_signal;

% =========================================================================
% DEFINE THE ULTRASOUND TRANSDUCER
% =========================================================================
% physical properties of transducer (in m)
element_width = 0.470e-3;
element_length = 14e-3;
kerf_width = 0.007e-3;

% physical properties of the transducer (in grid points)
transducer.number_elements = 64;    % total number of transducer elements
transducer.element_width = round(element_width/dx)    % width of each element 
transducer.element_length = round(element_length/dx)  % length of each element
transducer.element_spacing = round(kerf_width/dx)     % spacing (kerf width) between the elements
transducer.radius = inf;            % radius of curvature of the transducer [m]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements*transducer.element_width ...
    + (transducer.number_elements - 1)*transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
transducer.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);

% properties used to derive the beamforming delays
transducer.sound_speed = 1540;              % sound speed [m/s]
transducer.focus_distance = 30e-3;          % focus distance [m]
transducer.elevation_focus_distance = 48.8e-3;% focus distance in the elevation plane [m]
transducer.steering_angle = 0;              % steering angle [degrees]

% apodization
transducer.transmit_apodization = 'Rectangular';    
transducer.receive_apodization = 'Rectangular';

% define the transducer elements that are currently active
transducer.active_elements = ones(transducer.number_elements, 1);

% append input signal used to drive the transducer
transducer.input_signal = input_signal;

% create the transducer using the defined settings
transducer = makeTransducer(kgrid, transducer);

% print out transducer properties
transducer.properties;

% =========================================================================
% DEFINE SENSOR MASK
% =========================================================================

% define a sensor mask through the central plane
sensor.mask = zeros(Nx, Ny, Nz);
switch MASK_PLANE
    case 'xy'
        % define mask
        sensor.mask(:, :, Nz/2) = 1;
        
        % store y axis properties        
        Nj = Ny;
        j_vec = kgrid.y_vec;
        j_label = 'y';
        
    case 'xz'
        % define mask
        sensor.mask(:, Ny/2, :) = 1;
        
        % store z axis properties
        Nj = Nz;
        j_vec = kgrid.z_vec;
        j_label = 'z';
        
    case 'xy-xz'
        sensor.mask(:, Ny/2, :) = 1;
        sensor.mask(:, :, Nz/2) = 1;

    case 'onaxis'
        sensor.mask(:, Ny/2, Nz/2) = 1;
    
    case 'vol'
        % only look at axial positions +/- 5 mm from the focal depth
        axial_bound = 5e-3; % mm
        axial_indices = round(Nx*((transducer.focus_distance-axial_bound)/x)):...
                        round(Nx*((transducer.focus_distance+axial_bound)/x));
                        
        % only look at middle 25% of lateral position data.
        lat_indices = round((3/8)*Ny):round((5/8)*Ny);
        
        % look at middle 50% of elevational position data.
        ele_indices = round(Nz/4):round(3*Nz/4);
        sensor.mask(axial_indices, lat_indices, ele_indices) = 1;
    
    % add quarter symmetric case here
end 


% set the record mode such that only the rms and peak values are stored
if USE_STATISTICS
    sensor.record = {'p_rms', 'p_max'};
else
    sensor.record = {'p', 'u_non_staggered'};
end

% =========================================================================
% RUN THE SIMULATION
% =========================================================================

% set the input settings
input_args = {'DisplayMask', transducer.all_elements_mask, ...
              'PMLInside', false, ...
              'PlotPML', false, ...
              'PMLSize', [PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE], ...
              'DataCast', DATA_CAST, ...
              'DataRecast', true, ...
              'PlotScale', [-source_strength/2, source_strength/2]
              };

% stream the data to disk in blocks of 100 if storing the complete time
% history 
%{
if ~USE_STATISTICS
    input_args = [input_args {'StreamToDisk', 100}];
end
%}

if ~exist(h5output, 'file'),
    % run the simulation
    %sensor_data = kspaceFirstOrder3DC(kgrid, medium, transducer, sensor, input_args{:});
    sensor_data = kspaceFirstOrder3D(kgrid, medium, transducer, sensor, input_args{:}, 'SaveToDisk', h5input);

else,
    disp('Reading in sensor data... ')
    sensor_data = h5read(h5output, '/p');
    disp('Done reading sensor data.')
    % =========================================================================
    % COMPUTE THE BEAM PATTERN USING SIMULATION STATISTICS
    % =========================================================================

    if USE_STATISTICS
        
        % reshape the returned rms and max fields to their original position
        sensor_data.p_rms = reshape(sensor_data.p_rms, [Nx, Nj]);
        sensor_data.p_max = reshape(sensor_data.p_max, [Nx, Nj]);
        
        % plot the beam pattern using the pressure maximum
        figure;
        imagesc(j_vec*1e3, (kgrid.x_vec - min(kgrid.x_vec(:)))*1e3, sensor_data.p_max/1e6);
        xlabel([j_label '-position [mm]']);
        ylabel('x-position [mm]');
        title('Total Beam Pattern Using Maximum Of Recorded Pressure');
        %colormap(jet(256));
        c = colorbar;
        ylabel(c, 'Pressure [MPa]');
        axis image;
        
        % plot the beam pattern using the pressure rms
        figure;
        imagesc(j_vec*1e3, (kgrid.x_vec - min(kgrid.x_vec(:)))*1e3, sensor_data.p_rms/1e6);
        xlabel([j_label '-position [mm]']);
        ylabel('x-position [mm]');
        title('Total Beam Pattern Using RMS Of Recorded Pressure');
        %colormap(jet(256));
        c = colorbar;
        ylabel(c, 'Pressure [MPa]');
        axis image;

        % end the example
        return
        
    end

    % =========================================================================
    % COMPUTE THE BEAM PATTERN FROM THE AMPLITUDE SPECTRUM
    % =========================================================================

    % reshape the sensor data to its original position so that it can be
    % indexed as sensor_data(x, j, t)
    sensor_data = reshape(sensor_data, [Nx, Nj, kgrid.Nt]);

    % compute the amplitude spectrum
    [freq, amp_spect] = spect(sensor_data, 1/kgrid.dt, 'Dim', 3);

    % compute the index at which the source frequency and its harmonics occur
    [f1_value, f1_index] = findClosest(freq, tone_burst_freq);
    [f2_value, f2_index] = findClosest(freq, tone_burst_freq*2);

    % extract the amplitude at the source frequency and store
    beam_pattern_f1 = amp_spect(:, :, f1_index);

    % extract the amplitude at the second harmonic and store
    beam_pattern_f2 = amp_spect(:, :, f2_index);       

    % extract the integral of the total amplitude spectrum
    beam_pattern_total = sum(amp_spect, 3);

    % plot the beam patterns
    figure;
    imagesc(j_vec*1e3, (kgrid.x_vec - min(kgrid.x_vec(:)))*1e3, beam_pattern_f1/1e6);
    xlabel([j_label '-position [mm]']);
    ylabel('x-position [mm]');
    title('Beam Pattern At Source Fundamental');
    colormap(jet(256));
    c = colorbar;
    ylabel(c, 'Pressure [MPa]');
    axis image;
    print('-dpng', 'beam_pattern_1harm.png');

    figure;
    imagesc(j_vec*1e3, (kgrid.x_vec - min(kgrid.x_vec(:)))*1e3, beam_pattern_f2/1e3);
    xlabel([j_label '-position [mm]']);
    ylabel('x-position [mm]');
    title('Beam Pattern At Second Harmonic');
    colormap(jet(256));
    c = colorbar;
    ylabel(c, 'Pressure [kPa]');
    axis image;
    print('-dpng', 'beam_pattern_2harm.png');

    figure;
    imagesc(j_vec*1e3, (kgrid.x_vec - min(kgrid.x_vec(:)))*1e3, beam_pattern_total/1e6);
    xlabel([j_label '-position [mm]']);
    ylabel('x-position [mm]');
    title('Total Beam Pattern Using Integral Of Recorded Pressure');
    colormap(jet(256));
    c = colorbar;
    ylabel(c, 'Pressure [MPa]'); 
    axis image;
    print('-dpng', 'beam_pattern.png');

    % =========================================================================
    % PLOT DIRECTIVITY PATTERN AT FOCUS
    % =========================================================================

    % compute the directivity at each of the harmonics
    directivity_f1 = squeeze(beam_pattern_f1(round(transducer.focus_distance/dx), :));
    directivity_f2 = squeeze(beam_pattern_f2(round(transducer.focus_distance/dx), :));

    % normalize
    directivity_f1 = directivity_f1./max(directivity_f1(:));
    directivity_f2 = directivity_f2./max(directivity_f2(:));

    % compute relative angles from transducer
    if strcmp(MASK_PLANE, 'xy')
        horz_axis = ((1:Ny) - Ny/2)*dy;
    else
        horz_axis = ((1:Nz) - Nz/2)*dz;
    end
    angles = 180*atan2(horz_axis, transducer.focus_distance)/pi;

    % plot the directivity
    figure;
    plot(angles, directivity_f1, 'k-', angles, directivity_f2, 'k--');
    axis tight;
    set(gca, 'FontSize', 12);
    xlabel('Angle [deg]');
    ylabel('Normalized Amplitude (dB)');
    legend('Fundamental', 'Second Harmonic', 'Location', 'NorthWest'); 
    print('-dpng', 'directivity.png');
end

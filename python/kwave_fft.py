from time import time

import numpy as np
import h5py
import scipy.io as sio
import matplotlib.pyplot as plt



def read_data(h5file):
    return h5py.File(h5file, 'r')


def get_freq_vector(mat_contents, kwave_data):
    '''
    Reads in number of time points (Nt), sample period (dt), and
    FFT frequency vector (based on Nyquist sampling limit) from the
    given .h5 file. Stores these variables in the mat_contents
    dictionary.
    
    INPUTS:
    mat_contents - dict() containing sim variables to be saved to
                   a .MAT file.
    kwave_data - h5py.File() containing data stored in the .h5 file.
    
    OUTPUTS:
    mat_contents - updated dict() with time and frequency variables.
    '''
    mat_contents['Nt'] = kwave_data['Nt'][0, 0, 0]
    mat_contents['dt'] = kwave_data['dt'][0, 0, 0] # s 
    mat_contents['fs'] = 1/mat_contents['dt'] # Hz
    
    # Frequency vector will be cut to only include up to the 2nd harmonic.
    f = np.linspace(0, mat_contents['fs'], mat_contents['Nt']) # Hz
    max_harmonic = 2
    max_freq = (max_harmonic+1.5)*mat_contents['f0'] # Hz
    # Find index that contains the maximum frequency we are interested in,
    # then shorten frequency vector to contain values below the max frequency.
    mat_contents['max_f_index'] = np.argmax(f>max_freq)
    mat_contents['f'] = f[:mat_contents['max_f_index']]
    
    return mat_contents


def map_intensity_to_planes(mat_contents, intensity, kwave_data):
    mat_contents['Nx'] = kwave_data['Nx'][0, 0, 0]
    mat_contents['Ny'] = kwave_data['Ny'][0, 0, 0]
    mat_contents['Nz'] = kwave_data['Nz'][0, 0, 0]

    intensity_lat_plane = np.zeros((mat_contents['Nx'],
                                    mat_contents['Ny']))
    intensity_ele_plane = np.zeros((mat_contents['Nx'],
                                    mat_contents['Nz']))              
                                                                   
    # Note: 0-based indexing vs. MATLAB's 1-based indexing.
    print 'Re-constructing lateral and elevational planes'
    indices = (mat_contents['Nx'], mat_contents['Ny'], mat_contents['Nz'])
    x, y, z = np.unravel_index(mat_contents['sensor_indices'], indices, order='F')
    for i in xrange(mat_contents['num_sensor_pts']):
        if z[i] == (mat_contents['Nz']/2)-1 and y[i] == (mat_contents['Ny']/2)-1:
            intensity_lat_plane[x, y] = intensity[i]
            intensity_ele_plane[x, z] = intensity[i]
        elif y[i] == (mat_contents['Ny']/2)-1:
            intensity_ele_plane[x, z] = intensity[i]
        elif z[i] == (mat_contents['Nz']/2)-1:
            intensity_lat_plane[x, y] = intensity[i]
        else:
            print 'ERROR: Cannot re-construct center planes.'
            break
        # Progress bar
        if (i % (mat_contents['num_sensor_pts']/100)) == 0:
            print '%.1f %% complete.' % (100*(i/mat_contents['num_sensor_pts']))
    return mat_contents
    
    
def main():
    mat_contents = dict()
    
    mat_contents['h5input'] = '/luscinia/nl91/kwave/scratch/pwr2p0/kwave_input_data.h5'
    mat_contents['h5output'] = '/luscinia/nl91/kwave/scratch/pwr2p0/kwave_output_data.h5'
    
    ## Open .h5 input and output file for reading.
    kwave_data = read_data(mat_contents['h5output'])
    kwave_input_data = read_data(mat_contents['h5input'])
    
    ## Set transducer excitation frequency here.
    #  - Need a better way to do this so we don't have to hard-code it each time
    mat_contents['f0'] = 2.36e6; # Hz
    
    ## Get frequency vector from number of sample points and sampling frequency.
    mat_contents = get_freq_vector(mat_contents, kwave_data)
    
    ## Get number of sensor points and sensor point indices.
    mat_contents['sensor_indices'] = kwave_input_data['sensor_mask_index'][0, 0, :]
    mat_contents['num_sensor_pts'] = mat_contents['sensor_indices'].size
    
    ## Calculate intensities by calculating the sum of squared pressures 
    ## over time.
    print 'Calculating intensities...'
    intensity = np.zeros((mat_contents['num_sensor_pts']))
    for i in xrange(mat_contents['Nt']):
        p_timestep = kwave_data['p'][0, i, :]
        intensity += np.square(p_timestep)
        if (i % (mat_contents['Nt']/20)) == 0:
            print '%.1f %% complete.' % (100*(i/mat_contents['Nt']))
    print ''
    
    ## Map intensity data to elevational and lateral planes.
    mat_contents = map_intensity_to_planes(mat_contents, intensity, kwave_data)

    ## Write data to .MAT file.
    sio.savemat('kwave_intensity.mat', mat_contents, oned_as='column')


if __name__ == '__main__':
    start = time()
    main()
    end = time()
    print 'Total runtime: %.2f ms' % ((end-start)*1000.0)

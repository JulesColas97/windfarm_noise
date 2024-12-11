import numpy as np
import h5py
import os

def delta_L_ip(path, file_name, freq, nx, nz):
    os.chdir(path)
    delta_L_freq = np.zeros((nx, nz, freq.shape[0]))
    print(delta_L_freq.shape)
    for temp in range(freq.shape[0]):
        delta_L_freq[:, :, temp] = np.transpose(np.array(h5py.File(file_name, 'r+')['solution']["%04d" % freq[temp]]['deltaL']))

    mesh = {}
    mesh['x_coord'] = np.transpose(np.array(h5py.File(file_name, 'r+')['mesh']['x']))
    mesh['z_coord'] = np.transpose(np.array(h5py.File(file_name, 'r+')['mesh']['y']))

    return delta_L_freq, mesh

def delta_L_ip_2dpe(path, file_name, freq, MX, MZ):
    os.chdir(path)
    delta_L_freq = np.zeros((MZ, MX, freq.shape[1]))
    for temp in range(freq.shape[1]):
        delta_L_freq[:, :, temp] = np.array(h5py.File(file_name, 'r+')["freq_%04d" % freq[0, temp]]).T
    return delta_L_freq

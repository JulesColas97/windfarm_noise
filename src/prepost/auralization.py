import numpy as np
import matplotlib.pyplot as plt
import os
from .spl import SplField
from scipy.fft import rfft, rfftfreq, irfft
from scipy.signal.windows import tukey
from scipy.interpolate import interp1d
from scipy.io.wavfile import write, read
from . import audio_tools

SAMPLING = 44100
LMAX = 10
N_SPH = (LMAX + 1)**2


def cart2sph(x, y, z, degree=True):
    x, y, z = np.array(x), np.array(y), np.array(z)
    r, theta, phi = np.zeros(x.shape), np.zeros(x.shape), np.zeros(x.shape)
    r = np.sqrt(x**2 + y**2 + z**2)  # Radius
    theta = np.arctan2(np.sqrt(x**2 + y**2), z)  # Azimuth
    phi = np.arctan2(y, x)  # Colatitude
    return r, theta, phi


def decode_room(signal: np.ndarray):
    print('decoding room ...')
    L = 3
    signal_listening_room = audio_tools.decoder(signal, lmax=L)
    return signal_listening_room


def decode_bin(signal: np.ndarray):
    L = 3
    signal_bin = audio_tools.decoder_binaural(signal, lmax=L)
    return signal_bin


def decode(signal: np.ndarray,lmax=3):
    print('decoding ...')
    N = signal.shape[1]
    fs = SAMPLING
    signal_listening_room = audio_tools.decoder(signal, lmax=lmax)
    # Get speakers configuration
    spkr = audio_tools.lmfa_xyz_spkr()
    _, th_spkr, ph_spkr = audio_tools.cart2sph(
        spkr[:, 0], spkr[:, 1], spkr[:, 2])

    # Initialize binaural stereo sound file
    signal_binaural = np.zeros((N, 2))
    # Sum the binaural signals from each loudspeaker
    for i_spkr in range(spkr.shape[0]):
        # Retrieve IRs from sofa file
        ir_l, ir_r = audio_tools.get_hrir(th_spkr[i_spkr], ph_spkr[i_spkr])
        # Convolve noise with IRs
        signal_binaural += audio_tools.binaural(
            signal_listening_room[i_spkr, :], ir_l, ir_r)

    # Normalization
    signal_binaural /= signal_binaural.max()

    return signal_listening_room, signal_binaural


def interpolate_frequency(f0: np.ndarray, s0: np.ndarray, nsamples: int, sampling=44100):
    dt = 1/sampling
    f1 = rfftfreq(nsamples, dt)
    interpoland = interp1d(f0, s0, 'linear', bounds_error=False, fill_value=0)
    s1 = interpoland(f1)
    return f1, s1


def compute_rec_angles(splField: SplField):
    x = splField.x_grid[:, 0, 0].reshape(-1, 1, 1, 1, 1, 1)
    y = splField.y_grid[0, :, 0].reshape(1, -1, 1, 1, 1, 1)
    z = splField.z_grid[0, 0, :].reshape(1, 1, -1, 1, 1, 1)

    splField.wt.computeBeta()
    beta = splField.wt.beta.reshape(
        1, 1, 1, 1, splField.wt.Nblade, splField.wt.Nbeta)

    if splField.FULL_ROTATION:
        beta = np.concatenate((beta, beta+2*np.pi/3, beta+4*np.pi/3), axis=5)

    seg = splField.wt.seg.reshape(1, 1, 1, splField.wt.Nseg, 1, 1)

    # compute segment location
    zS = (np.cos(beta) * seg + splField.wt.href)
    xS = -np.sin(beta) * seg * np.sin(splField.wt.tau) + splField.xS
    yS = -np.sin(beta) * seg * np.cos(splField.wt.tau) + splField.yS
    # compute distance between source and receiver

    r, theta, phi = cart2sph((x-xS), (y-yS), (zS-z))
    return r, theta, phi


def compute_receiver_time(splField: SplField):
    print('compute receiver time ...')
    c0 = 343
    # (nx,ny,nz,src.wt.Nseg,src.wt.Nblade,len(freq),src.Nbeta)
    x = splField.x_grid[:, 0, 0].reshape(-1, 1, 1, 1, 1, 1)
    y = splField.y_grid[0, :, 0].reshape(1, -1, 1, 1, 1, 1)
    z = splField.z_grid[0, 0, :].reshape(1, 1, -1, 1, 1, 1)

    splField.wt.computeBeta()
    beta = splField.wt.beta.reshape(
        1, 1, 1, 1, splField.wt.Nblade, splField.wt.Nbeta)
    beta = np.concatenate((beta, beta+2*np.pi/3, beta+4*np.pi/3), axis=5)
    seg = splField.wt.seg.reshape(1, 1, 1, splField.wt.Nseg, 1, 1)

    # compute segment location
    # modif source :
    zS = (np.cos(beta) * seg + splField.wt.href)
    xS = -np.sin(beta) * seg * np.sin(splField.wt.tau) + splField.xS
    yS = -np.sin(beta) * seg * np.cos(splField.wt.tau) + splField.yS

    # compute distance between source and receiver
    R = np.sqrt((x - xS)**2 + (y - yS)**2 + (z - zS)**2)

    # absolute time of signal reception
    splField.t = R/c0 + beta[:, :, :, :, 0,
                             :].reshape(1, 1, 1, 1, 1, -1)/splField.wt.omega

    # copy first angle at the end
    splField.t = np.concatenate(
        (splField.t, splField.t[..., 0:1]+2*np.pi/splField.wt.omega), 5)
    print('done.')


def convert_freq_to_time(f0: np.ndarray, s0: np.ndarray, nsamples: int):
    frequencies, spp = interpolate_frequency(f0, s0, nsamples)
    phase = np.random.rand(len(frequencies))*2*np.pi
    spectrum = np.sqrt(spp)*np.exp(1j*phase)
    # print(spectrum)
    signal = irfft(spectrum)
    return signal


def freq_to_time(splField: SplField, x: float, y: float, z: float,
                 sampling=44100, lmax=3):
    Nsph = (lmax + 1)**2
    ix = np.argmin(np.abs(splField.x_grid[:, 0, 0]-x))
    iy = np.argmin(np.abs(splField.y_grid[0, :, 0]-y))
    iz = np.argmin(np.abs(splField.z_grid[0, 0, :]-z))

    # compute angle between source and receiver
    r, theta_rad, phi_rad = compute_rec_angles(splField)
    splField.compute_real_receiver_time(loop=False, last=True)

    print('compute time signal  ...')
    if splField.third:
        print('signal in third octave not possible to process.')
        quit()
    spp = 10**(splField.SPL_seg[ix, iy, iz, :, :, :, :]/10)*2e-5

    r = r[ix, iy, iz, :, :, :]
    theta = theta_rad[ix, iy, iz, :, :, :] * 180 / np.pi
    phi = phi_rad[ix, iy, 0, :, :, :] * 180 / np.pi
    t = splField.t[ix, iy, iz, :, :, :]

    t = t - np.min(t)
    T = 2*np.pi/(splField.wt.omega)

    nsamples_tot = int(T*sampling)

    temp_signal = np.zeros((nsamples_tot))
    temp_signal_sph = np.zeros((Nsph, nsamples_tot))

    it = 0
    overlap = 20
    alpha = 0.1
    for iblade in range(spp.shape[1]):
        # for iblade in range(1):
        for iseg in range(spp.shape[0]):
            it = int(t[iseg, iblade, 0]*sampling)
            for ibeta in range(0, spp.shape[3]):
                deltaT = (t[iseg, iblade, ibeta+1]-t[iseg, iblade, ibeta])
                nsamples = int(sampling*deltaT) + 2*overlap
                i0 = it - overlap
                iend = it + overlap + int(deltaT*sampling)
                indices = np.arange(i0, iend, 1) % nsamples_tot
                # print(spp[iseg,iblade,:,ibeta])
                si = convert_freq_to_time(splField.frequencies,
                                          spp[iseg, iblade, :, ibeta],
                                          nsamples+1)

                si_sph = audio_tools.encoder(si, lmax,
                                             theta[iseg, iblade, ibeta],
                                             phi[iseg, iblade, ibeta])

                window = tukey(nsamples, alpha)
                temp_signal[indices] = si[:nsamples] * \
                    window + temp_signal[indices]
                temp_signal_sph[:, indices] = si_sph[:, :nsamples] * \
                    window[None, :] + temp_signal_sph[:, indices]
                it = it+int(deltaT*sampling)

    # temp_signal = 0.8*temp_signal/np.max(temp_signal)
    # temp_signal_sph = 0.8*temp_signal_sph/np.max(temp_signal_sph)
    print(temp_signal_sph.shape)

    return temp_signal_sph


def freq_to_spherical(splField: SplField, x: float, y: float, z: float, lmax: float = 3, sampling=44100):
    Nsph = (lmax + 1)**2
    ix = np.argmin(np.abs(splField.x_grid[:, 0, 0]-x))
    iy = np.argmin(np.abs(splField.y_grid[0, :, 0]-y))
    iz = np.argmin(np.abs(splField.z_grid[0, 0, :]-z))

    # compute angle between source and receiver

    r, theta_rad, phi_rad = compute_rec_angles(splField)
    splField.compute_real_receiver_time(loop=True, last=True)
    # splField.compute_real_receiver_time(loop=False, last=True)

    print('compute time signal  ...')
    fc = [50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000]
    Nfc = [1,  1,  1,   1,   1,  1,   2,   2,   3,   4,   4,   4,   5,  5]
    splField.compute_third_octave(fc, Nfc)

    spp = np.sum(10**(splField.SPL_seg[ix, iy, iz, :, :, :, :]/10)*4e-10,
                 axis=2)

    # spp = spp
    r = r[ix, iy, iz, :, :, :]
    theta = theta_rad[ix, iy, iz, :, :, :] * 180 / np.pi
    phi = phi_rad[ix, iy, 0, :, :, :] * 180 / np.pi
    t = splField.t[ix, iy, iz, :, :, :]

    print(theta.shape)
    # plt.plot(theta[0, 1, :])
    # plt.plot(theta[0, 2, :])
    # plt.figure()
    # plt.plot(phi[0, 0, :])
    # plt.plot(phi[0,1,:])
    # plt.plot(phi[0,2,:])
    # plt.show()
    # quit()

    T = 2*np.pi/(splField.wt.omega)
    nsamples_tot = int(T*sampling)
    temp_signal_sph = np.zeros((Nsph, nsamples_tot))

    it = 0
    overlap = 0
    alpha = 0.1
    print("test")
    print(spp.shape[1])
    for iblade in range(spp.shape[1]):
        # for iblade in range(1):
        for iseg in range(spp.shape[0]):
            t0 = t[iseg, iblade, 0]
            i_first = int((t0 % T) * sampling)
            it = i_first
            # it = int(t[iseg, iblade, 0]*sampling)
            for ibeta in range(0, spp.shape[2]):
                tend = t[iseg, iblade, ibeta+1]
                iend = int((tend % T) * sampling)
                if iend > it:
                    indices = np.arange(it, iend, 1)
                else:
                    indices = np.concatenate([np.arange(it, nsamples_tot, 1),
                                               np.arange(0, iend, 1)], axis=0)

                # deltaT = (t[iseg, iblade, ibeta+1]-t[iseg, iblade, ibeta])
                # nsamples = int(sampling*deltaT) + 2*overlap
                # i0 = it - overlap
                # iend = it + overlap + int(deltaT*sampling)
                # indices = np.arange(i0, iend, 1) % nsamples_tot

                spl_sph = audio_tools.encoder(spp[iseg, iblade, ibeta],
                                              lmax, theta[iseg, iblade, ibeta],
                                              phi[iseg, iblade, ibeta])
                # window = tukey(nsamples, alpha)
                # if overlap > 0:
                #     temp_signal_sph[:, indices] = spl_sph[:, :] * \
                #         window[None, :] + temp_signal_sph[:, indices]
                # else:
                #    temp_signal_sph[:, indices] = (spl_sph[:, :] +
                #                                   temp_signal_sph[:, indices])
                temp_signal_sph[:, indices] = (spl_sph[:, :] +
                                               temp_signal_sph[:, indices])
                t0 = tend
                it = iend
                # it = it+int(deltaT*sampling)

    # temp_signam_sph = 10*np.log10(temp_signal_sph/2e-5)

    # temp_signal = 0.8*temp_signal/np.max(temp_signal)
    # temp_signal_sph = 0.8*temp_signal_sph/np.max(temp_signal_sph)

    return temp_signal_sph


def normalized_signals(directory: str = './'):
    max_pressure = 0
    for file in os.listdir(directory):
        if file.startswith("bin"):
            rate, signal_bin = read(directory+file)
            signal_max = np.max(signal_bin)

            if signal_max > max_pressure:
                max_pressure = signal_max

    for file in os.listdir(directory):
        if file.startswith("bin"):
            rate, signal_bin = read(directory+file)
            signal_bin = signal_bin/max_pressure
            write(directory+file, rate, signal_bin)
    print(max_pressure)
    max_pressure = 0
    for file in os.listdir(directory):
        if file.startswith("multi"):
            rate, signal_multi = read(directory+file)
            signal_max = np.max(signal_multi)
            if signal_max > max_pressure:
                max_pressure = signal_max
    print(max_pressure)
    # quit()
    for file in os.listdir(directory):
        if file.startswith("multi"):
            rate, signal_multi = read(directory+file)
            signal_multi = signal_multi/max_pressure
            write(directory+file, rate, signal_multi)


def combine_turbines(directory, length_s=10):

    coordinates = []
    turbines = []
    for file in os.listdir(directory):
        data = file.split('_')
        coordinates.append((int(data[2]), int(data[3][:-4])))
        turbines.append((data[1]))
    coordinates = list(dict.fromkeys(coordinates))
    turbines = list(dict.fromkeys(turbines))
    if 'sum' in turbines:
        turbines.remove('sum')
    print(turbines)

    nsamples = 44100*length_s
    for (x, y) in coordinates:
        signal = np.zeros((nsamples, 16))
        for t in turbines:
            print('multi_%s_%s_%s.wav' % (t, x, y))
            rate, signal_i = read(directory + 'multi_%s_%s_%s.wav' % (t, x, y))
            nsamples_i = signal_i.shape[0]
            print(nsamples_i)
            n = nsamples//nsamples_i
            r = nsamples % nsamples_i
            for ii in range(n):
                signal[ii*nsamples_i:(ii+1)*nsamples_i, :] += signal_i
            signal[n*nsamples_i:, :] += signal_i[:r, :]

        write(directory + 'multi_sum_%s_%s.wav' % (x, y), rate, signal)


def decode_room(signal: np.ndarray):
    print('decoding room ...')
    L = 3
    signal_listening_room = audio_tools.decoder(signal, lmax=L)
    return signal_listening_room


def decode_bin_multi(directory):
    L = 3
    for file in os.listdir(directory):
        if file.startswith("multi"):
            rate, signal = read(directory + file)
            print(signal.shape)
            signal_bin = audio_tools.decoder_binaural(signal.T, lmax=L)
            print(signal_bin.shape)
            write(directory + 'bin' + file[5:], rate, signal_bin)

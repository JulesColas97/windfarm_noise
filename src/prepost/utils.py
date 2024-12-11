import numpy as np
# import tikzplotlib
import os
import scipy.spatial.qhull as qhull
import itertools
import logging
import matplotlib.pyplot as plt


def uneven_tile(A: np.ndarray, reps: tuple):

    d = len(reps)
    if (d < A.ndim):
        tup = (1,)*(A.ndim-d) + reps

    even_reps = tuple(int(rep) for rep in reps)
    rest_reps = tuple(int((rep-int(rep))*n) for rep, n in zip(reps, A.shape))
    shape_out = tuple((s*et + rt)
                      for s, et, rt in zip(A.shape, even_reps, rest_reps))
    c = np.tile(A, even_reps)

    if c.size > 0:
        for idim, rep in enumerate(rest_reps):
            if rep > 0:
                indices = np.arange(0, rep)
                c = np.concatenate(
                    (c, np.take(c, indices, axis=idim)), axis=idim)

    return c.reshape(shape_out)


def uneven_loop(A: np.ndarray, N: int):
    n = A.shape[-1]
    if N == n:
        return A
    if N < n:
        print("cant loop because new Nt smaller than old Nt")
        return -1
    ncopy = N // n
    nrest = N % n
    new_shape = list(A.shape)
    new_shape[-1] = N
    new_shape = tuple(new_shape)

    A_new = np.zeros(new_shape)

    for ii in range(ncopy):
        A_new[..., ii*n:(ii+1)*n] = A
    if nrest > 0:
        A_new[..., (ncopy)*n:] = A[..., :nrest]
    return A_new


def computeThirdOctaveFrequencies(fc, Nfc):
    if len(fc) != len(Nfc):
        print('Fc and Nfc must be of same length')
        return
    Nfreq = sum(Nfc)
    print(Nfreq)
    freq = np.array([])
    for ii in range(len(fc)):
        freq = np.concatenate(
            (freq, np.round(fc[ii]*2**((2*np.arange(1, Nfc[ii]+1)/(Nfc[ii]+1) - 1)/6))))
    print(len(freq))
    return (freq)


def Aweight(f):
    Af = 12200**2*f**4./(f**2+20.6**2)/(f**2+12200**2) / \
        (f**2+107.7**2)**0.5/(f**2+737.9**2)**0.5
    dBA = 20*np.log10(Af/0.7943)
    return dBA

# def chkList(lst):
#     return len(set(lst)) == 1


def chkList(lst: list) -> (bool, int):
    ele = lst[0]
    # Comparing each element with first item
    for ii, item in enumerate(lst):
        if not np.all(ele == item):
            return False, ii
    return True, 0


def interp_weights(xy, uv, d=2):
    tri = qhull.Delaunay(xy)
    simplex = tri.find_simplex(uv)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uv - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))


def interpolate(values, vtx, wts):
    return np.einsum('nj,nj->n', np.take(values, vtx), wts)


def atm_absorption(T0, p0, rh, f):
    # function return the atmospheric attenuation of sound
    # due to the thermo-viscous effects and relaxation of oxygen and nitrogen.
    #
    # Usage: [alpha] = atm_absorption(T0,p0,rh,f)
    #         alpha - attenuation of sound for input parameters in dB/m
    #         T0 - temperature in K
    #         p0 - static pressure in pascal
    #         rh - relative humidity en #
    #         f - frequency of sound (may be a vector)
    #
    #
    # References:   Salomons p.109-111

    p0_ref = 1.01325e+05  # reference static pressure (pa)

    T0_tpw = 273.15  # triple point in K
    T0_ref = 293.15  # ref temp in K

    rho = p0/p0_ref
    tau = T0/T0_ref

    # calculate saturation pressure
    Csat = -6.8346*(T0_tpw/T0)**1.261 + 4.6151
    p0_sat = p0_ref*10**Csat
    h = rh*p0_sat/p0  # absolute humidity

    # Scaled relaxation frequency for Nitrogen
    frN = rho*tau**(-1/2)*(9 + 280*h*np.exp(-4.17*(tau**(-1/3)-1)))

    # scaled relaxation frequency for Oxygen
    frO = rho*(24 + 40400*h*(0.02+h)/(0.391+h))

    # attenuation coefficient in dB/m
    b1 = 0.1068*np.exp(-3352/T0)/(frN + f**2/frN)
    b2 = 0.01275*np.exp(-2239.1/T0)/(frO + f**2/frO)
    alpha = 8.686*f**2*tau**(1/2)*(1.84e-11/rho + tau**(-3)*(b1 + b2))
    return alpha


def cos_window(n, overlap):
    w = np.ones((n,))
    if overlap == 0:
        return w
    if overlap > n/2:
        logging.warning('error overlap larger than signal')
        quit()
    w[:overlap] = np.sin(0.5*np.arange(0, overlap, 1)*np.pi/(overlap-1))**2
    w[-overlap:] = np.cos(0.5*np.arange(0, overlap, 1)*np.pi/(overlap-1))**2

    return w

# Integration of third octave band spectrum s of size (...,Nf)
# with len(frequencies)=Nf


def integrateThirdOctave(frequencies, s):
    # #name of third octave bands
    # one_third_freq_preferred = np.array([5, 6.3, 8, 10, 12.5, 16, 20, 25,  31.5,
    # 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000])
    one_third_freq_preferred = np.array(
        [50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000])

    # check if frequencies and spectrum are broadcastable
    if (len(frequencies) != s.shape[-1]):
        print('last dimension of s must be equal to number of frequencies')
        return (one_third_freq_preferred, -1)
    # number of band
    Nband = len(one_third_freq_preferred)
    # df = frequencies[1] - frequencies[0]
    # real central frequency of each band
    one_third_bands = np.zeros((Nband, 2))
    one_third_freq = 1000*((2**(1/3)))**(np.arange(1, Nband+1)-Nband)

    # upper and lower limit of each band
    one_third_bands[:, 0] = one_third_freq/(2**(1/6))
    one_third_bands[:, 1] = one_third_freq*(2**(1/6))
    # allocate memory for third octave band spevctrum
    bands = np.zeros(s.shape[0:-1] + (Nband,))
    # loop over the thir octave bands
    for ii in range(Nband):
        # find frequencies index inside the band
        idx = np.squeeze(np.array(np.nonzero((frequencies >= one_third_bands[ii, 0]) & (
            frequencies < one_third_bands[ii, 1]))))

        # if band outside frequency range
        if idx.size == 0:
            bands[..., ii] = 0
        # if only one frequency inside the band
        elif (idx.size == 1):
            df = one_third_freq[ii]*0.232
            bands[..., ii] = df*(s[..., idx])
        # I need to take a look at what I have done here (from matlab code)
        else:
            df = one_third_freq[ii]*0.232/len(idx)
            if np.amin(idx) == 0:
                bands[..., ii] = df*(0.75*s[..., idx[1]] +
                                     np.sum(s[..., idx[1:]], -1) +
                                     0.75*s[..., idx[-1]] + 0.25*(s[..., idx[-1]+1]))

            elif np.amax(idx) == s.shape[-1]-1:
                bands[..., ii] = df*(0.25*s[..., idx[0]-1] + 0.75*s[..., idx[0]]
                                     + np.sum(s[..., idx[1:-1]], -1)
                                     + 0.75*s[..., idx[-1]])
            else:
                bands[..., ii] = df*(0.25*s[..., idx[0]-1] + 0.75*s[..., idx[0]]
                                     + np.sum(s[..., idx[2:-1]], -1)
                                     + 0.75*s[..., idx[-1]] + 0.25*s[..., idx[-1]+1])

    return (one_third_freq_preferred, bands)


def save_figure(filepath: str = None, fig='gcf', clean: bool = True):
    plt.savefig(filepath) 
#     from matplotlib.lines import Line2D
#     from matplotlib.legend import Legend
#     Line2D._us_dashSeq = property(lambda self: self._dash_pattern[1])
#     Line2D._us_dashOffset = property(lambda self: self._dash_pattern[0])
#     Legend._ncol = property(lambda self: self._ncols)
#
#     if clean:
#         tikzplotlib.clean_figure()
#
#     isExist = os.path.exists(os.path.dirname(filepath))
#     if not isExist:
#         os.makedirs(os.path.dirname(filepath))
#     tikzplotlib.save(filepath=filepath, figure=fig,
#                      wrap=False,
#                      add_axis_environment=True,
#                      tex_relative_path_to_data='\\figpath',
#                      dpi=300,
#                      externalize_tables=True,
#                      override_externals=True,
#                      extra_axis_parameters={
#                          # 'axis background/.style={fill=gray!60}',
#                          'xtick distance= \\xtick',
#                          'ytick distance = \\ytick',
#                          'colorbar=\\colorb',
#                          'colorbar style ={title=\\clegend,at={(1.02,1)},anchor=north west,width=0.15cm,major tick length=0.15cm,tick pos=right}',
#                          'height=\\height',
#                          'width=\\width',
#                          'tick align=outside',
#                          'xmin=\\xmin',
#                          'xmax=\\xmax',
#                          'ymin=\\ymin', 'ymax=\\ymax',
#                          'xlabel=\\xlabel',
#                          'ylabel =\\ylabel'}
#                      )
#     print(tikzplotlib.Flavors.latex.preamble())
#
def save_simple_figure(filepath: str = None, fig='gcf', clean: bool = True):
    plt.savefig(filepath) 
#     from matplotlib.lines import Line2D
#     from matplotlib.legend import Legend
#     Line2D._us_dashSeq = property(lambda self: self._dash_pattern[1])
#     Line2D._us_dashOffset = property(lambda self: self._dash_pattern[0])
#     Legend._ncol = property(lambda self: self._ncols)
#     if clean:
#         tikzplotlib.clean_figure()
#
#     isExist = os.path.exists(os.path.dirname(filepath))
#     if not isExist:
#         os.makedirs(os.path.dirname(filepath))
#     tikzplotlib.save(filepath=filepath, figure=fig,
#                      wrap=False,
#                      add_axis_environment=True,
#                      tex_relative_path_to_data='\\figpath',
#                      dpi=300,
#                      externalize_tables=True,
#                      override_externals=True,
#                      extra_axis_parameters={
#                          # 'axis background/.style={fill=gray!60}',
#                          'height=\\height',
#                          'width=\\width',
#                          "axis on top",
#                          'colorbar style ={title=\\clegend,at={(1.02,1)},anchor=north west,width=0.15cm,major tick length=0.15cm,tick pos=right}',
#                          }
#                      )
#     print(tikzplotlib.Flavors.latex.preamble())

if __name__ == "__main__":

    window = cos_window(20, 10)

    itrans = [0, 13, 20, 25, 30]
    nsamples_tot = 30
    coef = [1, 2, 3, 4, 5]
    overlap = 3

    it = itrans[0]

    final = np.zeros((30,))

    for ii in range(len(itrans)-1):
        deltai = itrans[ii+1] - itrans[ii]
        i0 = it - overlap
        iend = it + overlap + deltai
        indices = np.arange(i0, iend, 1) % nsamples_tot
        final[indices] = coef[it, None]

    plt.plot(final)
    plt.show()

    quit()

    test0 = np.zeros((30,))
    test1 = np.zeros((30,))
    test2 = np.zeros((30,))

    test2[:20] += window**2
    test0[:20] += window**2
    test2[10:] += window**2
    test1[10:] += window**2

    plt.plot(window)
    print(window)

    plt.figure()
    plt.plot(test0)
    plt.plot(test1)
    plt.plot(test2)
    plt.show()

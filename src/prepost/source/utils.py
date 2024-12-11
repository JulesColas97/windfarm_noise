from scipy.interpolate import CubicSpline
from scipy.interpolate import RegularGridInterpolator
import numpy as np

# interpolate 2D profile over blade segements position


def interp_atmos_data(z, u, znew):
    cs = CubicSpline(z, u)
    return cs(znew)


def interp_3D_atmos_data(atmos, xnew, ynew, znew):
    # check:
    if xnew.shape != ynew.shape or xnew.shape != znew.shape:
        raise Exception("new points must have same shape")
    points = np.concatenate(
        (znew.reshape(-1, 1), ynew.reshape(-1, 1), xnew.reshape(-1, 1)), axis=1)

    interpoland = RegularGridInterpolator(
        (atmos.z, atmos.y, atmos.x), atmos.u, method='linear',bounds_error=False)
    u_new = interpoland(points)
    u_new = u_new.reshape(xnew.shape)
    interpoland = RegularGridInterpolator(
        (atmos.z, atmos.y, atmos.x), atmos.epsilon, method='linear',bounds_error=False)
    epsilon_new = interpoland(points)
    epsilon_new = epsilon_new.reshape(xnew.shape)
    return u_new, epsilon_new

# go from (x,yz) coordinate system to (x,z,tau)


def cotte_conv(coord):
    receivers = {}
    receivers['x_rec'] = np.sqrt(coord['x_coord'] ** 2 + coord['y_coord'] ** 2).reshape(
        coord['x_coord'].shape[0] * coord['x_coord'].shape[1], 1)
    receivers['z_rec'] = (coord['z_coord']).reshape(
        coord['x_coord'].shape[0] * coord['x_coord'].shape[1], 1)
    receivers['tau'] = (np.arctan2(coord['y_coord'], coord['x_coord'])).reshape(
        coord['x_coord'].shape[0] * coord['x_coord'].shape[1], 1)
    return receivers


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


def c_round(temp):
    return round(temp * 100) / 100


def find_index(aoa, AOA):
    for temp in range(len(AOA)):
        if aoa <= AOA[temp]:
            break
    return temp

# for the definition of function which calculates the distance between each blade segment and the receiver locations


def R1_func(coord, H_ref):
    R1_array = np.sqrt((0 - coord['x_coord']) **
                       2 + (H_ref - coord['z_coord']) ** 2)
    return R1_array

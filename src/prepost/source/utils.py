from scipy.interpolate import CubicSpline
from scipy.interpolate import RegularGridInterpolator
import numpy as np

# interpolate 2D profile over blade segements position


def interp_atmos_data(z: np.ndarray, u: np.ndarray, znew: np.ndarray) -> np.ndarray:
    """
    Interpolate 1D atmospheric data using cubic splines.

    Args:
        z (np.ndarray): Original altitude values.
        u (np.ndarray): Original wind speed values.
        znew (np.ndarray): New altitude values for interpolation.

    Returns:
        np.ndarray: Interpolated wind speed values at new altitudes.
    """
    cs = CubicSpline(z, u)
    return cs(znew)


def interp_3D_atmos_data(atmos, xnew: np.ndarray, ynew: np.ndarray, znew: np.ndarray) -> tuple:
    """
    Interpolate 3D atmospheric data for given spatial coordinates.

    Args:
        atmos (Atmosphere): Atmosphere object containing wind speed and turbulence data.
        xnew (np.ndarray): New x-coordinates.
        ynew (np.ndarray): New y-coordinates.
        znew (np.ndarray): New z-coordinates.

    Returns:
        tuple: Interpolated wind speed and turbulence dissipation values.
    """
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


def cotte_conv(coord: dict) -> dict:
    """
    Convert (x, y, z) coordinate system to (x, z, tau) system.

    Args:
        coord (dict): Dictionary containing 'x_coord', 'y_coord', and 'z_coord'.

    Returns:
        dict: Converted coordinates with 'x_rec', 'z_rec', and 'tau'.
    """
    receivers = {}
    receivers['x_rec'] = np.sqrt(coord['x_coord'] ** 2 + coord['y_coord'] ** 2).reshape(
        coord['x_coord'].shape[0] * coord['x_coord'].shape[1], 1)
    receivers['z_rec'] = (coord['z_coord']).reshape(
        coord['x_coord'].shape[0] * coord['x_coord'].shape[1], 1)
    receivers['tau'] = (np.arctan2(coord['y_coord'], coord['x_coord'])).reshape(
        coord['x_coord'].shape[0] * coord['x_coord'].shape[1], 1)
    return receivers


def computeThirdOctaveFrequencies(fc: list, Nfc: list) -> np.ndarray:
    """
    Compute third-octave band center frequencies.

    Args:
        fc (list): List of center frequencies.
        Nfc (list): List of frequency counts for each band.

    Returns:
        np.ndarray: Computed third-octave frequencies.
    """
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


def Aweight(f: np.ndarray) -> np.ndarray:
    """
    Compute A-weighting for given frequencies.

    Args:
        f (np.ndarray): Array of frequencies in Hz.

    Returns:
        np.ndarray: A-weighted values in decibels.
    """
    Af = 12200**2*f**4./(f**2+20.6**2)/(f**2+12200**2) / \
        (f**2+107.7**2)**0.5/(f**2+737.9**2)**0.5
    dBA = 20*np.log10(Af/0.7943)
    return dBA


def c_round(temp: float) -> float:
    """
    Round a number to two decimal places.

    Args:
        temp (float): Number to be rounded.

    Returns:
        float: Rounded value to two decimal places.
    """
    return round(temp * 100) / 100


def find_index(aoa: float, AOA: np.ndarray) -> int:
    """
    Find the index of the given angle of attack (aoa) in the array AOA.

    Args:
        aoa (float): Angle of attack to search for.
        AOA (np.ndarray): Array of available angles.

    Returns:
        int: Index of the found angle in the AOA array.
    """
    for temp in range(len(AOA)):
        if aoa <= AOA[temp]:
            break
    return temp

def R1_func(coord: dict, H_ref: float) -> np.ndarray:
    """
    Calculate the distance between each blade segment and the receiver locations.

    Args:
        coord (dict): Dictionary containing 'x_coord' and 'z_coord'.
        H_ref (float): Reference height for the source.

    Returns:
        np.ndarray: Array of computed distances.
    """
    R1_array = np.sqrt((0 - coord['x_coord']) **
                       2 + (H_ref - coord['z_coord']) ** 2)
    return R1_array

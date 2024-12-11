#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 7 2023
@author: Loïc Alexandre

---------------------------------
Python tools for spatial audio
---------------------------------

"""

# Packages -------------------------------------------------------------------------------------------

import numpy as np
import sofa
# import sounddevice as sd
import pyshtools
from os import path
from scipy.signal import fftconvolve
# Utils -------------------------------------------------------------------------------------------

def c(celerity=343):

    """Return the celerity of sound in m/s.
         
    Returns
    -------
    c : float
        The sound celerity in m/s.
    
    Parameters
    ----------
    celerity : float, default=343
        The celerity in m/s.
    """

    c = celerity
    return c

def k(f, c=343):
    
    """Return the acoustic wave number in rad/m.
         
    Returns
    -------
    k : float
        The acoustic wave number in rad/m.
    
    Parameters
    ----------
    f : float
        The frequency in Hz.
    c : float
        The sound speed in m/s (default = 343 m/s).
    """
    return 2 * np.pi * f / c

def omega(f):

    """Return the acoustic pulsating frequency in rad.
         
    Returns
    -------
    omega : float
        The acoustic pulsating frequency in rad.
    
    Parameters
    ----------
    f : float
        The frequency in Hz.
    """
    return 2 * np.pi * f

def cart2sph(x, y, z, degree=True):
    """Return the spherical coordinate ``r, theta, phi`` of an array ``x, y, z`` in 
    cartesian coordinates.
    
    Returns
    -------
    r : float
        The value in spherical coordinates
    theta : float
        The value in spherical coordinates
    phi : float
        The value in spherical coordinates
    
    Parameters
    ----------
    x : float
        The value in Cartesian coordinate
    y : float
        The value in Cartesian coordinate
    z : float
        The value in Cartesian coordinate
    degree : bool, default=True
        If true, returns the spherical coordinates in degree.
        
    Description
    -----------
    Spherical coordinate system used:        
        Define :math:`r` to be distance (radius) from a point to the origin,
        :math:`\\theta` to be the colatitude angle and  :math:`\\phi` to be the 
        azimuthal angle. The ISO 80000-2:2019 convention is used.
    
    """
    x, y, z = np.array(x), np.array(y), np.array(z)
    if x.shape != y.shape != z.shape:
        raise ValueError('Arg arrays must have the same shape')
    r, theta, phi = np.zeros(x.shape), np.zeros(x.shape), np.zeros(x.shape)
    r = np.sqrt(x**2 + y**2 + z**2) # Radius
    theta = np.arctan2(np.sqrt(x**2 + y**2), z) # Azimuth
    phi = np.arctan2(y, x) # Colatitude 
    
    
    if degree==True:
        return r, np.rad2deg(theta), np.rad2deg(phi)
    else:
        return r, theta, phi


def sph2cart(r, theta, phi, deg=True):
    
    """Return the Cartesian coordinate ``x, y, z`` of an array ``r, theta, phi`` in 
    spherical coordinates.
    
    Returns
    -------
    x : float
        The value in Cartesian coordinate
    y : float
        The value in Cartesian coordinate
    z : float
        The value in Cartesian coordinate
        
    Parameters
    ----------
    r : float
        The value in spherical coordinates
    theta : float
        The value in spherical coordinates
    phi : float
        The value in spherical coordinates
        
    Description
    -----------
    Spherical coordinate system used:        
        Define :math:`r` to be distance (radius) from a point to the origin,
        :math:`\\theta` to be the colatitude angle and  :math:`\\phi` to be the 
        azimuthal angle. The ISO 80000-2:2019 convention is used.
    
    """
    r, theta, phi = np.array(r), np.array(theta), np.array(phi)
    if (r.shape != theta.shape != phi.shape):
        raise ValueError('Arg arrays must have the same shape')
    x, y, z = np.zeros(r.shape), np.zeros(r.shape), np.zeros(r.shape)
    
    if deg==True:
        x = r * np.sin(np.deg2rad(theta)) * np.cos(np.deg2rad(phi))
        y = r * np.sin(np.deg2rad(theta)) * np.sin(np.deg2rad(phi))
        z = r * np.cos(np.deg2rad(theta))
    else:
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)

    return x, y, z

# LMFA -------------------------------------------------------------------------------------------

def lmfa_room_dim():

    """Returns the x, y and z dimensions of the listening room at LMFA.
            
    Returns
    -------
    room_dim : array (3,)
        The dimensionss of the listening room.
        
    """ 
    
    room_dim = np.array([3, 4, 2.7])
    
    return room_dim


def lmfa_xyz_spkr():
    
    """Returns the x, y and z coordinates of the listening room loudspeakers at LMFA.
            
    Returns
    -------
    xyz : array (22, 3)
        The coordiantes of the loudspeakers.
        
    """ 
    
    xyz = np.array([[1.4,1.9,0],[1.4,1,0],[1.4,0,0],[1.4,-1,0],[1.4,-1.9,0], # Center front
                        [0,1.9,0],[0,-1.9,0], # Center center
                        [-1.4,1,0],[-1.4,0,0],[-1.4,-1,0], # Center back
                        [1.4,1.9,1.25],[1.4,0,1.25],[1.4,-1.9,1.25], # Up front
                        [0,1.9,1.25],[0,0,1.25],[0,-1.9,1.25], # Up center
                        [-1.4,1,1.25],[-1.4,0,1.25],[-1.4,-1,1.25], # Up back
                        [1.4,1.9,-1.25],[1.4,0,-1.25],[1.4,-1.9,-1.25]]) # Bottom front
    
    return xyz

# Binaural -------------------------------------------------------------------------------------------

def get_hrir(theta, phi):

    """
    Return the Head Related Impulse Responses corresponding from the closest angles ``theta, phi`` as FIR filters.
    
    Returns
    -------
    rir_l : array
        The left ear RIR 
    rir_r : float
        The right ear RIR
        
    Parameters
    ----------
    theta : float
        The azimuthal angle in degrees
    phi : float
        The colatitude angle in degrees
        
    Description
    -----------
    Spherical coordinate system used:        
        Define :math:`r` to be distance (radius) from a point to the origin,
        :math:`\\theta` to be the azimuthal angle in the :math:`xy`-plane from 
        the :math:`x`-axis with :math:`0\leq\\theta<2\pi`, and :math:`\phi` to 
        be the polar angle (also known as the zenith angle and colatitude, with
        :math:`\phi=90^\circ - \delta` where :math:`\delta` is the latitude) 
        from the positive :math:`z`-axis with :math:`0\leq \phi \leq \pi`. 
        This is the convention commonly used in mathematics. 
    
    """
    HRTF_path = path.dirname(__file__) + '/HRTF_44100.sofa' # Path to sofa HRIR file
    HRTF = sofa.Database.open(HRTF_path) # Open sofa data

    # Retrieve source positions of the HRIR measurements
    src_positions = HRTF.Source.Position.get_values(system="spherical") 

    # Get cartesian coordinates of the source positions for : 
    # r=1, colatitude = 90 - latitude, azimuth
    vec_cart_src = np.stack(sph2cart(1, 90-src_positions[:,1], src_positions[:,0]), axis=1)
    # Get cartesian coordinates of the given source angles
    vec_cart_source = np.stack(sph2cart(1, theta, phi), axis=0)
    # Retrieve the index of the nearest grid point based on dot product
    idx = np.argmax(vec_cart_source@vec_cart_src.T).squeeze()
    #
    # n_samples = int(np.floor(((src_positions[0,2]-0.1)/c())*HRTF.Data.SamplingRate.get_values(indices={"M":idx})))
    # Get IR for the left ear (receiver 0)
    ir_l = HRTF.Data.IR.get_values(indices={"M":idx, "R":0, "E":0})
    # ir_l = np.concatenate((np.roll(ir_l, -n_samples)[:-n_samples], np.zeros(n_samples)))
    # Get IR for the right ear (receiver 1)
    ir_r = HRTF.Data.IR.get_values(indices={"M":idx, "R":1, "E":0})
    # ir_r = np.concatenate((np.roll(ir_r, -n_samples)[:-n_samples], np.zeros(n_samples)))

    return ir_l, ir_r


def binaural(sound, ir_l, ir_r):

    """
    Return normalized binaural stereo sound.
    
    Returns
    -------
    binaural_sound : array
        The binaural stereo sound
        
    Parameters
    ----------
    sound : array
        The sound file to convolve
    rir_l : array
        The IR filter for the left ear
    rir_r : array
        The IR filter for the right ear
    fs : float
        The sampling frequency
    play : bool, default=True
        If true, plays the binaural sound file using sounddevice.
        
    Description
    -----------
    
    """

    # Lenght of the sound file
    N = sound.shape[0]
    # Initialize binaural stereo sound matrix
    binaural_sound = np.zeros((N, 2))
    # Convolve left side
    binaural_sound[:,0] = fftconvolve(sound, ir_l)[:N]
    # Convolve right side
    binaural_sound[:,1] = fftconvolve(sound, ir_r)[:N]

    return binaural_sound


def decoder_binaural(blm: np.ndarray, lmax: int=3, method:str='mode_matching', src='plane'):

    """Return decoded sound file.

    Returns
    -------
    binaural_sound : 2-D array
        sound.shape[0] by 2 matrix 
        
    Parameters
    ----------
    blm : array
        The encoded sound file of shape lmax by N samples
    lmax : int
        The ambisonic order (lmax > 0)
    method : string, default='mode_matching'
        Deconding method. Should be 'mode_matching'
    src : string, default='plane'
        The type of source to consider for encoding. Should be 'plane or 'monopole'.
        
    Description
    -----------

    """

    HRTF_path = path.dirname(__file__) + '/HRTF_44100.sofa' # Path to sofa HRIR file
    HRTF = sofa.Database.open(HRTF_path) # Open sofa data

    # Retrieve source positions of the HRIR measurements
    src_positions = HRTF.Source.Position.get_values(system="spherical") 

    if method=='mode_matching':

        C = np.zeros(shape=((lmax + 1)**2, src_positions.shape[0]))
        
        if src=='plane':
                   
            # Spherical harmonic coefficients of the loudspeakers
            for src_count in range(src_positions.shape[0]):       
                C[:, src_count] = spharm_acn(lmax, 90-src_positions[src_count,1], src_positions[src_count,0])

        elif src=='monopole' :

            print('Sorry, not implemented yet')

        elif (src != 'plane' != 'monopole'):
            raise ValueError('src must be plane or monopole')

        # Matrix pseudo inverse to obtain the decoding matrix
        decoder = np.linalg.pinv(C, rcond=1/100)
        # Calculate the loudspeaker signals from the target SH sound field
        w = decoder@blm
        # Initialize binaural stereo sound file
        binaural_sound = np.zeros((blm.shape[1], 2))
        # Sum the binaural signals from each loudspeaker
        for i_src in range(src_positions.shape[0]):
            # Retrieve IRs from sofa file
            ir_l, ir_r = get_hrir(90-src_positions[i_src,1], src_positions[i_src,0])
            # Convolve noise with IRs
            
            binaural_sound += binaural(w[i_src,:], ir_l, ir_r)

        return binaural_sound


# Spherical harmonics -------------------------------------------------------------------------------------------

def lm2acn(coeffs):
    """ 
    Out of Class version of lm2acn.
    Reshapes a matrix of coefficients from i-l-m to ACN 
    (Ambisonic Channel Number ordering).

    """
    lmax = coeffs.shape[2]-1
    newcoeffs = np.zeros([(lmax+1)**2], dtype = type(coeffs[0, 0, 0]))
    for acn in range((lmax+1)**2):
        l = int(np.sqrt(acn))
        m = int(acn - l**2 - l)
        if m<0:
            i = 1
        else:
            i = 0
        newcoeffs[acn] = coeffs[i, l, abs(m)]
    return newcoeffs


def acn2lm(coeffs):
    """ 
    Converts the matrix of coefficients of an SHCoeffs instance 
    from ACN (Ambisonic Channel Number ordering) to i-l-m
    """
    lmax = int(np.sqrt(coeffs.size))
    newcoeffs = np.zeros([2, lmax, lmax], dtype = type(coeffs[0]))
    for acn in range(coeffs.size):
        l = int(np.sqrt(acn))
        m = int(acn - l**2 - l)
        if m<0:
            i = 1
        else:
            i = 0
        newcoeffs[i, l, np.abs(m)] = coeffs[acn]
    coeffs = newcoeffs
    return coeffs

def spharm_acn(lmax, theta, phi, **kwargs):
    """
    Returns the :math:`4\pi`-normalized real spherical harmonic of degree ``l`` 
    and order ``m`` in the ACN ordering.

    Returns
    -------
    y : float
        The spherical harmonic ``ylm`` evaluated at argument (``theta``, 
        ``phi``).

    Parameters
    ----------
    l : integer
        The spherical harmonic order.
    theta : float
        The colatitude angle in degree.
    phi : float
        The azimuth angle in degree.
        
    Description
    -----------
    The :math:`4\pi`-normalized real spherical harmonic of degree l and order m
    is defiened as follows:

    .. math ::
        \begin{equation}
        Y_{lm}(\theta,\varphi) = \sqrt{(2l+1)(2 - \delta(m)) 
        \frac{(l-|m|)!}{(l+|m|)!}} 
        \left \lbrace \begin{array}{ll} P_{lm}
        (\sin(\varphi)) \cos (m \theta) & \mbox{if $m \ge 0$} \\
        P_{l|m|}(\sin (\varphi)) \sin (|m| \theta) & \mbox{if $m < 0$},
        \end{array} \right.
        \end{equation}

    where :math:`P_{lm}` is the associated Legendre polynomial of degree
    :math:`l` and order :math:`m`.
    """
    Ylm = pyshtools.expand.spharm(lmax, theta, phi, **kwargs)
    Yacn = lm2acn(Ylm)
    return Yacn

def encoder(sound, lmax, theta, phi, src='plane'):

    """Return encoded sound file.

    Returns
    -------
    blm : 2-D array
        (lmax+1)**2 by sound.shape[0] matrix
        
    Parameters
    ----------
    sound : array
        The sound file to convolve
    lmax : int
        The ambisonic order
    theta : float
        The colatitude angle of the source
    phi : float
        the azimuthal angle of the source
    src : string, default='plane'
        The type of source to consider for encoding. Should be 'plane or 'monopole'.
        
    Description
    -----------

    """
    # Get the spherical harmonic coefficients at the given angles
    ylm = spharm_acn(lmax, theta, phi)

    if src=='plane':
        # Encoding the source signal
        blm = np.tile(sound, ((lmax+1)**2,1)).T * ylm

        return blm.T

    elif src=='monopole':
        print('Sorry, not implemented yet ...')

    elif (src != 'plane' != 'monopole'):
        raise ValueError('src must be plane or monopole')

        
def decoder(blm, lmax, spkr=None, method='mode_matching', src='plane'):

    """Return decoded sound file.

    Returns
    -------
    w : 2-D array
        spkr.shape[0] by sound.shape[0] matrix
        
    Parameters
    ----------
    blm : array
        The encoded sound file
    lmax : int
        The ambisonic order (lmax > 0)
    spkr : array, default=None
        The N loudspeaker positions in cartesian coordinates (N,3). If None, the loudspeaker configuration from LMFA is used.
    method : string, default='mode_matching'
        Deconding method. Should be 'mode_matching'
    src : string, default='plane'
        The type of source to consider for encoding. Should be 'plane or 'monopole'.
        
    Description
    -----------

    """

    if spkr is None:
        spkr = lmfa_xyz_spkr()
    
    # Get spherical coordinates of loudspeakers (no need for r values)
    _, th_spkr, ph_spkr = cart2sph(spkr[:,0], spkr[:,1], spkr[:,2])

    if method=='mode_matching':

        C = np.zeros(shape=((lmax + 1)**2, spkr.shape[0]))
        
        if src=='plane':
                   
            # Spherical harmonic coefficients of the loudspeakers
            for spkr_count in range(spkr.shape[0]):       
                C[:, spkr_count] = spharm_acn(lmax, th_spkr[spkr_count], ph_spkr[spkr_count])
            
            # Matrix pseudo inverse to obtain the decoding matrix
            d = np.linalg.pinv(C, rcond=1/100)
            # Calculate the loudspeaker signals from the target SH sound field
            w = d@blm

            return w

        elif src=='monopole' :

            print('Sorry, not implemented yet')

        elif (src != 'plane' != 'monopole'):
            raise ValueError('src must be plane or monopole')

import numpy as np

def make_gaussian(x, w0):
    # Purely transverse 1D scalar field
    return np.exp(-(x**2 / w0**2)) + 0j

def calc_sag(r, R, k):
    if R == 0: 
        return np.zeros_like(r)
    c = 1.0 / R
    return (c * r**2) / (1.0 + np.sqrt(1.0 - (1.0 + k) * c**2 * r**2))

def make_lens_slice(x, z_cur, R1, k1, R2, k2, thick, diam, n_lens, n_bg):
    r = np.abs(x)
    z_front = calc_sag(r, R1, k1)
    z_back = thick + calc_sag(r, R2, k2)
    
    in_lens = (r <= diam/2.0) & (z_cur >= z_front) & (z_cur <= z_back)
    
    n_map = np.full_like(r, n_bg, dtype=float)
    n_map[in_lens] = n_lens
    return n_map

def wpm_step(E_in, kx, k0, n_slice, dz, window):
    E_out = np.zeros_like(E_in, dtype=complex)
    E_fft = np.fft.fft(E_in)
    k_trans2 = kx**2

    for nm in np.unique(n_slice):
        mask = (n_slice == nm)
        kz2 = (k0 * nm)**2 - k_trans2
        kz2[kz2 < 0] = 0.0 
        kz = np.sqrt(kz2)
        
        # Angular spectrum of plane waves propagation
        E_prop_fft = E_fft * np.exp(1j * kz * dz)
        E_prop = np.fft.ifft(E_prop_fft)
            
        E_out += mask * E_prop
        
    return E_out * window

def calc_overlap_1d(E_sim, x, w_out):
    E_target = np.exp(-(x / w_out)**2) + 0j
    num = np.abs(np.sum(np.conj(E_sim) * E_target))**2
    den = np.sum(np.abs(E_sim)**2) * np.sum(np.abs(E_target)**2)
    if den == 0: return 0.0
    
    eta_1d = num / den
    # Square the 1D overlap to approximate the 2D cylindrical fiber coupling efficiency
    eta_2d = eta_1d**2 
    return eta_2d
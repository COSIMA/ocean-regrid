import numpy as np

#pythran export apply_weights(float64[:,:], (int, int), int, int, int32[:], int32[:], float64[:])
def apply_weights(src, dest_shape, n_s, n_b, row, col, s):
    """
    Apply ESMF regridding weights.
    """

    dest = np.zeros(dest_shape).flatten()
    src = src.flatten()

    weighted_src =  s[0:n_s]*src[col[0:n_s]-1]

    for i in range(0, n_s):
        dest[row[i]-1] = dest[row[i]-1] + weighted_src[i]

    return dest.reshape(dest_shape)

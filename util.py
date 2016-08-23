
import numpy as np

def normalise_lons(lons, data=None):
    """
    Normalise longitudes to 0-360 deg. Perform the same transformation on data.
    """

    # Remove -ves
    new_lons = np.copy(lons)
    new_lons[lons < 0] = lons[lons < 0] + 360
    lons_copy = np.copy(new_lons)

    if data is not None:
        new_data = np.copy(data)
    else:
        new_data = None

    # Use changes in 2nd derivative to find jumps. Then offset (the +2) to get
    # element directly after jump. It just works OK.
    jumps = list(np.where(np.diff(new_lons[0, :], 2))[0][::2] + 2)

    # Beginning and sizes of continuous segments of lons.
    segs = [0] + jumps
    sizes = np.diff(segs + [len(new_lons[0, :])])

    # Sort according to value of lon at segment begin index.
    segs = zip(new_lons[0, segs], segs, sizes)
    src_segs = sorted(segs, key=lambda x : x[0])

    dest_idx = 0
    for i, (_, src_idx, size) in enumerate(src_segs):
        new_lons[:, dest_idx:dest_idx+size] = lons_copy[:, src_idx:src_idx+size]

        if new_data is not None:
            if len(data.shape) == 3:
                new_data[:, :, dest_idx:dest_idx+size] = data[:, :, src_idx:src_idx+size]
            else:
                new_data[:, dest_idx:dest_idx+size] = data[:, src_idx:src_idx+size]

        dest_idx += size

        if i+1 == len(segs):
            break

    return new_lons, new_data

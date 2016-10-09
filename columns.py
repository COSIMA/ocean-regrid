
def regrid_columns(src_data, src_grid, dest_grid, temp_or_salt):
    """
    Regrid vertical columns of data from src to dest.

    1. Fill in source bathymetry with nearest neighbours.
    2. interpolate and extrapolate down so that dest columns are filled to the
    bottom.
    """

    assert temp_or_salt == 'temp' or temp_or_salt == 'salt'
    assert len(src_data.shape) == 3
    assert src_data.shape[0] == len(src_grid.z)

    # This gets modified.
    src_data = np.ma.copy(src_data)
    src_data_a = np.ma.copy(src_data)
    src_data_b = np.ma.copy(src_data)

    # Create masked array of the correct shape.
    tmp = np.zeros((len(dest_grid.z), src_data.shape[1], src_data.shape[2]))
    new_data = np.ma.array(tmp, mask=np.ones_like(tmp), copy=True)

    # We need to fill in missing values in the source bathymetry. Do this in
    # steps:
    #   1a. fill with horizontal nearest neighbour so that deep point
    #   have valued based on neighbours at the same depth.
    #   1b. fill with nearest neighbour from above in the same column. Before
    #   doing this make sure the top level is covered.
    #   1c. for each grid cell choose the 'best' one from the above with the
    #   intention of finding the most stable arrangement of the water column.

    # 1a. At every level fill everything with nearest neighbour. This
    # effectively removes bathymetry with the new (previously masked) deep
    # points having values based on neighbours at the same depth.
    for level in range(src_data.shape[0]):
        ind = nd.distance_transform_edt(src_grid.mask[level, :, :],
                                        return_distances=False,
                                        return_indices=True)
        tmp = src_data_b[level, :, :]
        tmp = tmp[tuple(ind)]
        src_data_b[level, :, :] = tmp[:]

    # 1b. First fill in the top level.
    ind = nd.distance_transform_edt(src_grid.mask[0, :, :],
                                   return_distances=False,
                                   return_indices=True)
    tmp = src_data_a[0, :, :]
    tmp = tmp[tuple(ind)]
    src_data_a[0, :, :] = tmp[:]

    # Now fill in all missing values down the columns.
    for lat in range(src_data.shape[1]):
        for lon in range(src_data.shape[2]):
            ind = nd.distance_transform_edt(src_grid.mask[:, lat, lon],
                                            return_distances=False,
                                            return_indices=True)
        tmp = src_data_a[:, lat, lon]
        tmp = tmp[tuple(ind)]
        src_data_a[:, lat, lon] = tmp[:]

    # 1c. combine the best of above. I'm not sure what's 'best'. I guess deeper
    # water should be colder and more salty to make a more stable column.
    if temp_or_salt == 'temp':
        src_data[np.where(src_data_a[:] <= src_data_b[:])] = \
            src_data_a[np.where(src_data_a[:] <= src_data_b[:])]
        src_data[np.where(src_data_b[:] < src_data_a[:])] = \
            src_data_b[np.where(src_data_b[:] < src_data_a[:])]
    else:
        src_data[np.where(src_data_a[:] >= src_data_b[:])] = \
            src_data_a[np.where(src_data_a[:] >= src_data_b[:])]
        src_data[np.where(src_data_b[:] > src_data_a[:])] = \
            src_data_b[np.where(src_data_b[:] > src_data_a[:])]

    # Step 2. Iterate through columns and regrid each to the bottom of the destination column.
    for lat in range(src_data.shape[1]):
        for lon in range(src_data.shape[2]):
            # 1d linear interpolation/extrapolation
            new_data[:, lat, lon] = interp(dest_grid.z, src_grid.z, src_data[:, lat, lon])

    return new_data


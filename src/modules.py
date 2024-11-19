import numpy as np

def slice_coords(coords, dz):

    coords = np.asarray(coords)

    z_min = coords[:, 2].min()
    z_max = coords[:, 2].max()

    z_edges = np.arange(z_min, z_max + dz, dz)

    slice_indices = np.digitize(coords[:, 2], z_edges) - 1

    unique_indices = np.unique(slice_indices)

    slices = {
            idx: coords[slice_indices == idx]
            for idx in unique_indices
    }
    
    return slices, z_edges

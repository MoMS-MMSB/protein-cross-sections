import numpy as np
from tqdm import tqdm

def slice_coords(coords, n_slices):

    coords = np.asarray(coords)

    z_order = np.argsort(coords[:,2])
    sorted_coords = coords[z_order]
    # Split into n_slices groups of approximately equal size
    splits = np.array_split(sorted_coords, n_slices)
    
    # Create dictionary of slices
    slices = {i: split for i, split in enumerate(splits)}
    
    return slices

def results_by_slice(coords, n_slices):

    slices = slice_coords(coords, n_slices)
    results = {}
    for i in np.arange(len(slices)):
        centroid_xy = slices[i][:, :2].mean(axis=0)
        xy_distances =  np.sqrt(np.sum((slices[i][:, :2] - centroid_xy) ** 2, axis=1))

        results[i] = {
            'points': slices[i],
            'centroid': centroid_xy,
            'distances': xy_distances,
            'z_range': (slices[i][:, 2].min(), slices[i][:, 2].max)
            }
    return results

def distances_by_slice(coords, n_slices):

    slices = slice_coords(coords, n_slices)
    distances = []
    for i in np.arange(len(slices)):
        centroid_xy = slices[i][:, :2].mean(axis=0)
        xy_distances =  np.sqrt(np.sum((slices[i][:, :2] - centroid_xy) ** 2, axis=1))

        distances.append(np.mean(xy_distances))
    return distances

def distances_per_frame(u, n_slices):
    slice_distances = []
    for ts in tqdm(u.trajectory):
        protein = u.select_atoms("protein").positions
        slice_distances.append(distances_by_slice(protein, n_slices))

    avg_slice_distances = np.mean(slice_distances, axis=0)
    slice_std = np.std(slice_distances, axis=0)

    return avg_slice_distances, slice_std

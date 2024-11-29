import numpy as np
from tqdm import tqdm
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def slice_coords(coords, n_slices):

    coords = np.asarray(coords)

    z_order = np.argsort(coords[:,2])
    sorted_coords = coords[z_order]
    
    #print(int(np.max(coords[:, 2]) - np.min(coords[:, 2])))

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

## distances
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

## areas
def areas_by_slice(coords, n_slices):
    slices = slice_coords(coords, n_slices)
    areas = []
    for i in np.arange(len(slices)):
        try:
            hull = ConvexHull(slices[i][:, :2])
            areas.append(hull.area)
        except:
            areas.append(0)
    return areas 
 
def areas_per_frame(u, n_slices):
    slice_areas = []
    protein_lengths = []
    for ts in tqdm(u.trajectory):
        protein = u.select_atoms("protein").positions
        slice_areas.append(areas_by_slice(protein, n_slices))
        protein_lengths.append(np.max(protein[:, 2]) - np.min(protein[:, 2]))
    avg_slice_areas = np.mean(slice_areas, axis=0) * 10
    slice_std = np.std(slice_areas, axis=0) * 10
    avg_protein_length = np.mean(protein_lengths)
    return avg_slice_areas, slice_std, avg_protein_length

def length_per_frame(u):
    lengths = []
    for ts in tqdm(u.trajectory):
        protein = u.select_atoms("protein").positions
        lengths.append(np.max(protein[:, 2]) - np.min(protein[:, 2]))
    return np.average(lengths), np.std(lengths)

def plot_slice_hulls(coords, n_slices):
    slices = slice_coords(coords, n_slices)

    fig = plt.figure(figsize=(30,10))
    ax1 = fig.add_subplot(1, 3, 1)
    ax2 = fig.add_subplot(1, 3, 2,projection='3d')
    ax3 = fig.add_subplot(1, 3, 3)
    #ax1.view_init(0, -90, 0)
    ax2.view_init(45, 45, 0)
    #ax3.view_init(90, -90, 0)
    
    for i in np.arange(len(slices)):
        hull = ConvexHull(slices[i][:, :2])
        vertices = hull.points[hull.vertices]
        vertices_3d = np.column_stack([vertices, np.full(len(vertices), np.mean(slices[i][:, 2]))])
        vertices_3d = np.vstack([vertices_3d, vertices_3d[0]])
        
        ax1.plot(vertices_3d[:, 0], vertices_3d[:, 2],
            c=plt.cm.cool(i/n_slices), linewidth=2)
        ax2.plot(vertices_3d[:, 0], vertices_3d[:, 1], vertices_3d[:, 2],
            c=plt.cm.cool(i/n_slices))
        ax3.plot(vertices_3d[:, 0], vertices_3d[:, 1],
            c=plt.cm.cool(i/n_slices), linewidth=2)
    
    ax1.set_xlabel('X')
    ax1.set_ylabel('Z')
    
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')
    
    ax2.xaxis.pane.fill = False
    ax2.yaxis.pane.fill = False
    ax2.zaxis.pane.fill = False
    
    ax3.set_xlabel('X')
    ax3.set_ylabel('Y')

    plt.savefig("hull.png", bbox_inches="tight")

def voxelize_protein(protein_coords, box_dims, voxel_size=3):
    box_dims = box_dims[:3]
    n_voxels = np.ceil(np.array(box_dims[:3]) / voxel_size).astype(int)

    voxel_grid = np.zeros(n_voxels)

    voxel_indices = (protein_coords / voxel_size).astype(int)

    voxel_grid[voxel_indices[:, 0], voxel_indices[:, 1], voxel_indices[:, 2]] = 1

    occupied_area = np.sum(voxel_grid, axis=(0,1)) * (voxel_size ** 2)
    return occupied_area, voxel_grid
 
def voxels_per_frame(u, voxel_size):
    voxel_areas = []
    for ts in tqdm(u.trajectory):
        protein = u.select_atoms("protein").positions
        voxel_areas.append(voxelize_protein(protein, u.dimensions, voxel_size)[0])
    avg_voxel_areas = np.mean(voxel_areas, axis=0)
    voxel_std = np.std(voxel_areas, axis=0)
    return avg_voxel_areas, voxel_std


def plot_voxels(voxel_grid):
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')

    x, y, z = np.where(voxel_grid == 1)

    ax.voxels(voxel_grid, facecolors='blue', alpha=0.1, edgecolor='black')
    ax.set_xlabel('X (voxels)')
    ax.set_ylabel('Y (voxels)')
    ax.set_zlabel('Z (voxels)')
    ax.set_title('3D Voxel Representation of Protein')
    plt.savefig("voxels.png", bbox_inches="tight")

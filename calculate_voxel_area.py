import time
import modules
import argparse
import numpy as np
from tqdm import tqdm
import MDAnalysis as mda
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

parser.add_argument("-s", "--structure", help = "<.gro>  structure file", required = True)
parser.add_argument("-t", "--trajectory", help = "<.xtc>/<.trr> trajectory file")
parser.add_argument("-v", "--voxel_size", type=float, help = "voxel size in angstrom", default=3)
parser.add_argument("-o", "--output", default="voxel_area")
parser.add_argument("--plot_voxel", type=bool, help="return a nice 3d plot of the occupied voxels", default=False)
args = parser.parse_args()

try:
    if not args.trajectory:
        u = mda.Universe(args.structure)
    else:
        u = mda.Universe(args.structure, args.trajectory)
except NameError:
    print("Error loading the structure and/or trajectory file into MDAnalysis")

if args.trajectory:
    avg_slice_areas, slice_std, protein_length = modules.areas_per_frame(u, args.n_slices)
    protein = u.select_atoms("protein").positions
    z_range = np.arange(0, protein_length, protein_length/args.n_slices)
    print(len(avg_slice_areas), len(z_range))
    plt.errorbar(x=avg_slice_areas, y=z_range, xerr=slice_std,
                 c = "#d60036")
    plt.xlim(0, np.max(avg_slice_areas) + 50)
    results = np.stack([z_range, avg_slice_areas, slice_std]).T
    np.savetxt(f'{args.output}.csv', results, delimiter=',')
else:
    protein = u.select_atoms("protein").positions
    voxel_areas, voxel_grid = modules.voxelize_protein(protein, u.dimensions, args.voxel_size)
    plt.plot(voxel_areas, np.arange(0, len(voxel_areas)) * args.voxel_size)

plt.xlabel("Area (Å²)")
plt.ylabel("Z coords (Å)")
plt.title(f"{len(voxel_areas)} voxels in z- of {args.voxel_size} angstrom")
plt.savefig(f"{args.output}.png")

if args.plot_voxel:
    modules.plot_voxels(voxel_grid)

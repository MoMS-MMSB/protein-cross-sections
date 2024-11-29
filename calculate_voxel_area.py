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
    avg_voxel_areas, voxel_std = modules.voxels_per_frame(u, args.voxel_size)
    protein = u.select_atoms("protein").positions
    plt.errorbar(x=avg_voxel_areas, y=(np.arange(0, len(avg_voxel_areas))*3), xerr=voxel_std,
                 c = "#d60036")
    plt.xlim(-10, np.max(avg_voxel_areas) + 500)
    plt.title(f"{len(avg_voxel_areas)} voxels in z- of {args.voxel_size} angstrom")
    #results = np.stack([z_range, avg_slice_areas, slice_std]).T
    #np.savetxt(f'{args.output}.csv', results, delimiter=',')
else:
    protein = u.select_atoms("protein").positions
    voxel_areas, voxel_grid = modules.voxelize_protein(protein, u.dimensions, args.voxel_size)
    plt.plot(voxel_areas, np.arange(0, len(voxel_areas)) * args.voxel_size)
    plt.title(f"{len(voxel_areas)} voxels in z- of {args.voxel_size} angstrom")

plt.xlabel("Area (Å²)")
plt.ylabel("Z coords (Å)")
plt.savefig(f"{args.output}.png")

if args.plot_voxel:
    modules.plot_voxels(voxel_grid)

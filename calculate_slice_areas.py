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
parser.add_argument("-n", "--n_slices", type=int, help = "number of slices in z", default=50)
parser.add_argument("--plot_hull", type=bool, help="return a nice 3d plot of the slices and their convex hulls", default=False)

args = parser.parse_args()

try:
    if not args.trajectory:
        u = mda.Universe(args.structure)
    else:
        u = mda.Universe(args.structure, args.trajectory)
except NameError:
    print("Error loading the structure and/or trajectory file into MDAnalysis")

if args.trajectory:
    avg_slice_areas, slice_std = modules.areas_per_frame(u, args.n_slices)
    plt.errorbar(x=avg_slice_areas, y=np.arange(args.n_slices), xerr=slice_std,
                 c = "#d60036")
    plt.xlim(0, np.max(avg_slice_areas) + 50)
    results = np.stack([np.arange(args.n_slices), avg_slice_areas, slice_std]).T
    np.savetxt('test.csv', results, delimiter=',')
else:
    protein = u.select_atoms("protein").positions
    slice_areas = modules.areas_by_slice(protein, args.n_slices)
    plt.plot(slice_areas)

plt.savefig("test.png")

if args.plot_hull:
    modules.plot_slice_hulls(protein, args.n_slices)

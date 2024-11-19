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
parser.add_argument("-o", "--output", default="out")
args = parser.parse_args()

try:
    if not args.trajectory:
        u = mda.Universe(args.structure)
    else:
        u = mda.Universe(args.structure, args.trajectory)
except NameError:
    print("Error loading the structure and/or trajectory file into MDAnalysis")

if args.trajectory:
    avg_slice_distances, slice_std = modules.distances_per_frame(u, args.n_slices)
    plt.errorbar(x=avg_slice_distances, y=np.arange(args.n_slices), xerr=slice_std,
                 c = "#d60036")
    plt.xlim(0, np.max(avg_slice_distances) + 5)
    results = np.stack([np.arange(args.n_slices), avg_slice_distances, slice_std]).T
    np.savetxt(f"{args.output}.csv", results, delimiter=',')
else:
    protein = u.select_atoms("protein").positions
    slice_distances = modules.distances_by_slice(protein, args.n_slices)
    plt.plot(slice_distances)

plt.savefig(f"{args.output}.png")

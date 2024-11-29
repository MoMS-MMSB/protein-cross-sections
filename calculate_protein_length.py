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
parser.add_argument("-o", "--output", default="area")
args = parser.parse_args()

try:
    if not args.trajectory:
        u = mda.Universe(args.structure)
    else:
        u = mda.Universe(args.structure, args.trajectory)
except NameError:
    print("Error loading the structure and/or trajectory file into MDAnalysis")

if args.trajectory:
    length_avg, length_std = modules.length_per_frame(u)
    print(length_avg, length_std)
    #plt.errorbar(x=avg_slice_areas, y=z_range, xerr=slice_std,
    #             c = "#d60036")
    #plt.xlim(0, np.max(avg_slice_areas) + 50)
    #results = np.stack([z_range, avg_slice_areas, slice_std]).T
#    np.savetxt(f'{args.output}.csv', results, delimiter=',')
# else:
#     protein = u.select_atoms("protein").positions
#     slice_areas = modules.areas_by_slice(protein, args.n_slices)
#     plt.plot(slice_areas)
# plt.xlabel("Area (Å²)")
# plt.ylabel("Protein length in z (Å)")
# plt.title(f"Average area per slice, average std = {np.mean(slice_std):.2f}Å²")
# plt.savefig(f"{args.output}.png")
# 
# if args.plot_hull:
#     modules.plot_slice_hulls(protein, args.n_slices)

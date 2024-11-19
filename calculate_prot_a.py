import modules
import argparse
import MDAnalysis as mda

parser = argparse.ArgumentParser()

parser.add_argument("-s", "--structure", help = "<.gro>  structure file", required = True)
parser.add_argument("-t", "--trajectory", help = "<.xtc>/<.trr> trajectory file")
parser.add_argument("-dz", "--dz", help = "z interval for slicing (in Angstrom)", default=5)

args = parser.parse_args()

try:
    if not args.trajectory:
        u = mda.Universe(args.structure)
    else:
        u = mda.Universe(args.structure, args.trajectory)
except NameError:
    print("Error loading the structure and/or trajectory file into MDAnalysis")

protein = u.select_atoms("protein").positions
print(modules.slice_coords(protein, args.dz))

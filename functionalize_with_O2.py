# runs functionalize_catalyst.py, then add_o2_frozen.py
# creates a subdirectory of bare catalyst and another of O2-bound .xyz files
# see functionalize_catalyst.find_Hs() for list of accepted catalysts

import sys, os
import functionalize_catalyst, add_o2_frozen

catalyst = sys.argv[1]
core = sys.argv[2]
num_funcs = int(sys.argv[3])
num_molecules = int(sys.argv[4])

func_locs = functionalize_catalyst.find_Hs(catalyst)

outdir = functionalize_catalyst.run(catalyst, core, func_locs, num_funcs, num_molecules)
# catalyst, metal/core, expected number of funcs, number of molecules generated
# for nonmetal catalysts, you can put anything for the core (functionalize_catalyst.py ignores it)

for infile in os.listdir(outdir):
    molecule_O2 = add_o2_frozen.AddO2(outdir + "/" + infile)
    molecule_O2.run()

os.system("mkdir %s/catfunc/ %s/catfuncO2/" %(outdir, outdir))
os.system("mv %s/*O2*.xyz %s/catfuncO2/" %(outdir, outdir))
# won't work properly if catalyst name contains "O2"
os.system("mv %s/*.xyz %s/catfunc/" %(outdir, outdir))

# contains functions for randomly functionalizing atom indices
# ligand should be added to molSimplify library first (remove metal atom if applicable)
    # use the files on GitHub for non-porphyrin ligands
# there are currently no checks to make sure that number of requested items
    # does not exceed max possible number of items, so be careful (probably fine)

import os, sys
import molSimplify
from molSimplify.Scripts import structgen
from molSimplify.Classes import mol3D, atom3D
import numpy as np
import math, random

homedir = "/home/nricke/" # modify as needed
nonmetal_catalysts = ["mepyr", "tetrid", "tetry"]

def find_Hs(catalyst): # returns list of functionalizable atom indices
    # atom indices are one-indexed from .mol or .xyz file
    if catalyst == "porphyrin":
        return [9, 10, 17, 18, 24, 25, 30, 31, 33, 34, 35, 36]
    elif catalyst == "phencir":
        return [31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42]
    elif catalyst == "mepyr":
        return [19, 20, 25]
    elif catalyst == "tetrid":
        #return [30,31,34,35,36,37,38]
        return [31,32,35,36,37,38,39]
    elif catalyst == "tetry":
        return [16, 25, 26, 29, 30, 31, 32, 33]
    else:
        print("Not a valid ligand.")
        return None

def choose_items(list, num_items, no_repeats):
    chosen = [] # contains objects, not indices
    for n in range(num_items):
        index = random.randrange(0, len(list))
        if no_repeats: # list items should all be unique
            while list[index] in chosen:
                index = random.randrange(0, len(list))
        chosen.append(list[index])
    return chosen

def make_tempdir_outdir(catalyst, core):
    tempdir = homedir + "functionalizecatalysttempdir" # modify as needed
    while os.path.exists(tempdir):
        tempdir += "0"
    os.makedirs(tempdir)

    outdir = homedir + "%s%s_funcs" %(catalyst, core) # modify as needed; molSimplify expands ~ to home directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    else:
        version = 2
        outdir += str(version)
        while os.path.exists(outdir):
            version += 1
            outdir = homedir + "%s%s_funcs%d" %(catalyst, core, version)
        os.makedirs(outdir)

    return tempdir, outdir

def list_to_string(list): # list should not be empty
    if len(list) == 1:
        return str(list[0])
    rtn = "["
    for item in list:
        rtn += str(item) + "," # no spaces
    rtn = rtn.rstrip(",") + "]"
    return rtn

def functionalize(catalyst, core, possible_func_indices, expected_num_funcs, num_molecules, tempdir):
    # func_library: SMILES strings or built-in ligands; modify as needed
    # generally there are issues reading in strings containing []()
    func_library = ["B", "C", "N", "O", "F", "P", "S", "Cl", "Br", # basis set doesn't include I
                    "cyanide", "NC", "tricyanomethyl", "methylamine",
                    "dicyanamide", "nitroso",
                    "OC", "carboxyl",
                    "trifluoromethyl",
                    "thiocyanate",
                    "benzene_pi", "benzenethiol"]
    for n in range(1, num_molecules + 1):
        file_name = catalyst + core + "-functionalized" + str(n)
        if catalyst in nonmetal_catalysts:
            file_name += "TEMPMETAL"
        num_funcs = 0
        while (num_funcs < 1 or num_funcs > len(possible_func_indices)):
            num_funcs = np.random.poisson(expected_num_funcs)
        func_list = choose_items(func_library, num_funcs, False)
        func_indices = choose_items(possible_func_indices, num_funcs, True)

        ms_command = "molsimplify -core %s -oxstate 2 -coord 6 -geometry oct -spin 3 -lig %s -ligocc 1 -decoration %s -decoration_index %s -rundir %s -name %s" %(core, catalyst, list_to_string(func_list), list_to_string(func_indices), tempdir, file_name)
        print(ms_command)
        os.system(ms_command)

def collect_xyz_files(current_location, destination, delete_current_location):
    # janky workaround for molSimplify's directory structure
    # copies all .xyz files in current_location to destination
    # does not check for duplicate molecules or file names
    # should this skip "badjob" files? don't expect those to happen but just in case
    for subdir, dirs, files in os.walk(current_location):
        for file in files:
            if file.split('.')[-1] == "xyz":
                os.system("scp %s %s" %(os.path.join(subdir, file), destination))
    if delete_current_location:
        os.system("ls -R %s" %current_location)
        os.system("rm -R %s" %current_location)
        os.system("rm %sCLIinput.inp" %homedir)
        # CLIinput.inp actually gets saved to current working directory

def remove_TEMPMETAL(dir):
    dir = dir.rstrip("/")
    for file in os.listdir(dir):
        if file.split('.')[-1] != "xyz" or "TEMPMETAL" not in file:
            continue
        molecule = mol3D.mol3D()
        molecule.OBMol = molecule.getOBMol(dir + "/" + file, "xyz", ffclean = False)
        molecule.convert2mol3D()
        for i in molecule.findMetal(): # len(molecule.findMetal()) should initially be 1
            molecule.deleteatom(i)

        structgen.ffopt('UFF', molecule, [], 1, [], False, [], 200, False)
        # MMFF94/ffopt() in general can cause segfaults
        # you should still check the structures in Avogadro and optimize manually if necessary

        molecule.writexyz(dir + "/" + file[:file.index("TEMPMETAL")])
    # for file in os.listdir(dir):
    #     if "TEMPMETAL" in file:
    #         os.system("rm %s/%s" %(dir, file))

def run(catalyst, core, possible_func_indices, expected_num_funcs, num_molecules):
    if catalyst in nonmetal_catalysts:
        core = ""
        # for naming purposes (molSimplify will default to Fe, but it doesn't matter after remove_TEMPMETAL())
    tempdir, outdir = make_tempdir_outdir(catalyst, core)
    functionalize(catalyst, core, possible_func_indices, expected_num_funcs, num_molecules, tempdir)
    print("Functionalized .xyz files placed in " + outdir + "/")
    collect_xyz_files(tempdir, outdir, True)
    remove_TEMPMETAL(outdir)
    return outdir

if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2], find_Hs(sys.argv[1]), int(sys.argv[3]), int(sys.argv[4]))

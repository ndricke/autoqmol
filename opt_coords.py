# fix for optcdftsp bug in another script that ran cdftsp on unoptimized geometry
# this script takes optcdftsp.out file/directory and creates cdftsp.in using optimized geometry
# check file names first to prevent overwriting

import ChemData as CD
import os, sys

def run(infile): # optcdftsp.out
    # note that new .xyz and .in files will be placed in same directory as infile

    if infile.split('.')[-1] != "out":
        return None

    # write optimized coordinates to a new .xyz with the same file name
    xyz_fn = infile.split('.')[0] + ".xyz"
    xyz = CD.ConvergedCoordGrab(infile)
    if xyz != []:
        CD.WriteXYZ(xyz_fn, xyz)

#    # parse charge, mult, method, basis from infile
#    mol_flag = False
#    with open(infile, 'r') as f:
#        for line in f:
#            if "METHOD" in line:
#                method = line[6:].strip()
#            elif "BASIS" in line:
#                basis = line[5:].strip()
#            elif mol_flag:
#                mol_flag = False
#                try:
#                    charge, multiplicity = int(line.split(' ')[0]), int(line.split(' ')[-1])
#                except:
#                    continue
#            elif "$molecule" in line:
#                mol_flag = True

    # generate cdftsp.in # change qcMultIn.py path if necessary
    #os.system("python qcMultIn.py -f %s -c %d -m %d -method %s -basis %s -j cdftsp" %(xyz_fn, charge, multiplicity, method, basis))

if __name__ == "__main__":
    input = sys.argv[1].rstrip('/')
    if os.path.isfile(input):
        run(input)
    elif os.path.isdir(input):
        for file in os.listdir(input):
            run(input + '/' + file)
    else:
        print("Input should be a file or directory.")

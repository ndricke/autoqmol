# input: file or directory
# pickles Q-Chem .out files

import sys
import os
import Qdata
import pickle

def read_data(qchem_file): # does not check whether qchem_file is a .out
    qdata = Qdata.Qdata()
    qdata.readFile(qchem_file)
    return qdata

def pickle_data(in_filedir, pickle_fname):
    save_data = []
    for subdir, dirs, files in os.walk(in_filedir):
        for qfile in files:
            if qfile.split('.')[-1] == 'out':
                print(os.path.join(subdir,qfile))
                qdata = read_data(os.path.join(subdir,qfile))
                save_data.append(qdata)

    pickle.dump(save_data, open(pickle_fname, 'wb'))

if __name__ == "__main__":
    in_filedir = sys.argv[1]
    pickle_fname = sys.argv[2]
    pickle_data(in_filedir, pickle_fname)

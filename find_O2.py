# input: pickled catO2 directory or file
# identifies exactly two oxygen atoms as those in dioxygen; warning otherwise
# determines catalyst atom to which O2 is bound
# creates DataFrame containing data of interest for each catalyst, or
    # "NaN" if not applicable (the ints are atom indices)
# DataFrame indices now refer to list of atoms in Q-Chem .out file, not
    # atom_coords from pickle (Q-Chem index = atom_coords index + 1)

import numpy as np
import pandas as pd
import sys
import os
import Qdata
import pickle
from collections import OrderedDict

# used to create DataFrame
filenames, oxy1_list, oxy2_list, cat_list, energy_list = [], [], [], [], []
cat_O2_bond_lengths, oxy1_chelpgs, oxy2_chelpgs, total_oxy_chelpgs = [], [], [], []
cat_elements, cat_chelpgs = [], []

def original_index(O_index, atom_coords):
    # converts O_coords index to atom_coords index; precondition: O_index >= 0
    count = 0
    for i in range(len(atom_coords)):
        if atom_coords[i].split(' ')[0] == 'O':
            count += 1
        if count == O_index + 1:
            return i
    print(O_index, " is not a valid O_coords index.")

def get_coords(input_string):
    # input_string of the form "<atom> <x-coord> <y-coord> <z-coord>"
    return np.array([float(input_string.split(' ')[1]), float(input_string.split(' ')[2]), float(input_string.split(' ')[3])])

def distance(string1, string2):
    difference = get_coords(string1) - get_coords(string2)
    return abs(np.linalg.norm(difference))

def getKey(qdata):
    return qdata.filename

def order_loaded_pickle(inlist): # inlist should contain Qdatas
    # returns list entries, alphabetized by file name (starting from the last letter and going backwards)
    # for nonmetal catalysts, this will group spreadsheet lines by active site rather than functionalization
    # if you change my nonmetal O2-binding naming conventions, this might stop working
    inlist_ordered = inlist # inlist_ordered is returned; inlist is not modified
    for item in inlist_ordered:
        item.filename = item.filename[::-1]
    inlist_ordered = sorted(inlist_ordered, key = getKey)
    for item in inlist_ordered:
        item.filename = item.filename[::-1]
    return inlist_ordered

def create_df(in_pickle):
    save_data = order_loaded_pickle(pickle.load(open(in_pickle, "rb")))
    for qdata in save_data:
        try:
            atom_coords = qdata.listCoord()
        except:
            print("No coordinates in " + qdata.filename + ". Omitted from DataFrame.")
            continue

        filenames.append(qdata.filename)
        O2_found, multiple_O2 = False, False
        try:
            energy_list.append(qdata.E)
        except:
            energy_list.append(None)
        try:
            chelpgs = [float(item) for item in qdata.chelpg]
        except:
            chelpgs = []
            for i in atom_coords:
                chelpgs.append(None)

        O_coords = []
        for coord in atom_coords:
            if coord.split(' ')[0] == 'O':
                O_coords.append(coord)
        if len(O_coords) < 2:
            #print("Fewer than two oxygens in " + qdata.filename)
            oxy1_list.append(None)
            oxy2_list.append(None)
            oxy1_chelpgs.append(None)
            oxy2_chelpgs.append(None)
            total_oxy_chelpgs.append(None)
            cat_O2_bond_lengths.append(None)
            cat_list.append(None)
            cat_elements.append(None)
            cat_chelpgs.append(None)
            continue
        for index1 in range(len(O_coords) - 1):
            for index2 in range(index1 + 1, len(O_coords)):
                if distance(O_coords[index1], O_coords[index2]) < 1.6:
                    if O2_found == True:
                        multiple_O2 = True
                        break
                        # doesn't exit index1 for loop but still reduces unnecessary operations
                    else:
                        O2_found = True
                        oxy1, oxy2 = index1, index2 # O_coords index values for oxygens identified as O2
        if not O2_found:
            #print("No O2 detected in " + qdata.filename)
            oxy1_list.append(None)
            oxy2_list.append(None)
            oxy1_chelpgs.append(None)
            oxy2_chelpgs.append(None)
            total_oxy_chelpgs.append(None)
            cat_O2_bond_lengths.append(None)
            cat_list.append(None)
            cat_elements.append(None)
            cat_chelpgs.append(None)
            continue
        if multiple_O2:
            #print("Multiple O2 detected in " + qdata.filename)
            oxy1_list.append(None)
            oxy2_list.append(None)
            oxy1_chelpgs.append(None)
            oxy2_chelpgs.append(None)
            total_oxy_chelpgs.append(None)
            cat_O2_bond_lengths.append(None)
            cat_list.append(None)
            cat_elements.append(None)
            cat_chelpgs.append(None)
            continue

        shortest_dist, catalyst_index, oxygen_index = sys.float_info.max, sys.float_info.max, sys.float_info.max
        # catalyst_index and oxygen_index refer to atom_coords, not O_coords
        #oxygen_index unnecessary
        oxy1_orig, oxy2_orig = original_index(oxy1, atom_coords), original_index(oxy2, atom_coords)
        oxy1_list.append(oxy1_orig + 1)
        oxy2_list.append(oxy2_orig + 1)
        oxy1_chelpgs.append(chelpgs[oxy1_orig])
        oxy2_chelpgs.append(chelpgs[oxy2_orig])
        try:
            total_oxy_chelpg = chelpgs[oxy1_orig] + chelpgs[oxy2_orig]
        except:
            total_oxy_chelpg = None
        total_oxy_chelpgs.append(total_oxy_chelpg)
        # ensures oxy lists have correct length; last entry may be overwritten
        for atom_index in range(len(atom_coords)):
            if atom_index == oxy1_orig or atom_index == oxy2_orig:
                continue
            if distance(atom_coords[oxy1_orig], atom_coords[atom_index]) < shortest_dist:
                shortest_dist = distance(atom_coords[oxy1_orig], atom_coords[atom_index])
                catalyst_index = atom_index
                oxygen_index = oxy1_orig
                oxy1_list[-1] = oxy1_orig + 1
                oxy2_list[-1] = oxy2_orig + 1
                oxy1_chelpgs[-1] = chelpgs[oxy1_orig]
                oxy2_chelpgs[-1] = chelpgs[oxy2_orig]
                try:
                    total_oxy_chelpg = chelpgs[oxy1_orig] + chelpgs[oxy2_orig]
                except:
                    total_oxy_chelpg = None
                total_oxy_chelpgs[-1] = total_oxy_chelpg
            if distance(atom_coords[oxy2_orig], atom_coords[atom_index]) < shortest_dist:
                shortest_dist = distance(atom_coords[oxy2_orig], atom_coords[atom_index])
                catalyst_index = atom_index
                oxygen_index = oxy2_orig
                oxy1_list[-1] = oxy2_orig + 1
                oxy2_list[-1] = oxy1_orig + 1
                oxy1_chelpgs[-1] = chelpgs[oxy2_orig]
                oxy2_chelpgs[-1] = chelpgs[oxy1_orig]
                # oxy1 refers to the O bound to the catalyst, oxy2 the other O
                try:
                    total_oxy_chelpg = chelpgs[oxy1_orig] + chelpgs[oxy2_orig]
                except:
                    total_oxy_chelpg = None
                total_oxy_chelpgs[-1] = total_oxy_chelpg
        cat_O2_bond_lengths.append(shortest_dist)
        cat_list.append(catalyst_index + 1)
        cat_elements.append(atom_coords[catalyst_index].split(' ')[0])
        cat_chelpgs.append(chelpgs[catalyst_index])

    O2_dict = {'CatalystO2_File_Name': filenames, 'CatalystO2_Energy': energy_list, 'Active_Site': cat_list,
                                     'Active_Site_ID': cat_elements, 'Oxygen_1': oxy1_list, 'Oxygen_2': oxy2_list,
                                     'CatalystO2_Active_Site_CHELPG': cat_chelpgs, 'Oxygen_1_CHELPG': oxy1_chelpgs, 'Oxygen_2_CHELPG': oxy2_chelpgs,
                                     #'O2_CHELPG': total_oxy_chelpgs, # has rounding errors so do this manually on spreadsheet instead
                                     'Cat-O2_Bond_Length': cat_O2_bond_lengths}
    O2_df = pd.DataFrame.from_dict(OrderedDict(O2_dict))

    #print(O2_df)
    name = in_pickle.split('.')[0].split('/')[-1] # input pickle file name
    O2_df.to_csv(name + '.csv')
    O2_df.to_pickle(name + '_df.p')

if __name__ == "__main__":
    create_df(sys.argv[1])

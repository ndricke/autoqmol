#!/usr/bin/env python

"""
arguments: -
    
calls: ModificationDictionary, BindingSiteFragment
    
this script contains functions which pick what modifications to make on a mol object

!!! The weighted list has to be changed manually !!! (at the moment it is based on fg_summary.csv for 1-1006 sp
calculations) !!! (or it could be read directly from a fg_summary.csv)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from random import choice
import ModificationDictionary as md
from BindingSiteFragment import binding_site_fragment


weighted_list=['-F']*25+['-n']*18+['-OC']*5+['-NC2']*4+['-Br']*3+['-NH2']*3+['-NO2']+['-Cl']

def random_on_target(m,m_o2,cutoff,substruct=None):
    """
    arguments: mol object of cat, mol object of cat+o2, cutoff idx, substructure string (optional)
    
    returns: new mol object or None (if new mol object not possible)
        
    this function picks a random modification to be performed on the molecule for the given substructure 
    (if substruct is None, it picks a random substructure from the dictionary) on the ring next to the scaffold
    or on the binding site, neighbors or next neighbors
    """
    
    mw=Chem.RWMol(m)
    mw_o2=Chem.RWMol(m_o2)
    
    atom_list=sorted(list(binding_site_fragment(mw_o2)[-1]))
    atom_list=atom_list[:-2] #remove the two O atoms

    #Find the C neighbor that is part of the target ring (reactive fragment)
    for ngb in mw.GetAtomWithIdx(cutoff).GetNeighbors():
        #print(ngb.GetIdx(),ngb.GetSymbol())
        if ngb.GetSymbol()=='C' and ngb.GetIdx()>cutoff:
            ngb_idx=ngb.GetIdx()
            #print('the ngb idx is '+str(ngb_idx))
            break
    
    #Find the target ring
    for ring in mw.GetRingInfo().AtomRings():
        if (ngb_idx in ring) and (cutoff not in ring):
            target_ring=ring
            #print(target_ring)
            break
        
    atom_list.extend(target_ring)

    new_m=random_from_list(mw,atom_list,substruct=substruct)    
    
    return new_m

def general_modify(m, cutoff, substruct=None):

    mw=Chem.RWMol(m)
    atom_count = mw.GetNumAtoms()
    print("Atom count: ", atom_count)
    atom_list = range(cutoff, atom_count)
    print("Modif2 atom list: ", atom_list)

    new_m=random_from_list(mw, atom_list, substruct=substruct)    
    
    return new_m

def random_from_list(m,at_list,substruct=None,similarity=0):
    """
    arguments: mol object, list of atom indices, substructure string (optional), similarity index (optional)
    
    returns: new mol object or None (if new mol object not possible)
        
    this function picks a random modification to be performed on the molecule for the given substructure 
    (if substruct is None, it picks a random substructure from the dictionary) on a random atom from atom_list
    
    if similarity is set to a value between 0 and 1, it only returns modified structures that have a similarity 
    higher than or equal to that of the original structure (similarity evaluated using Morgan fingerprints)
    """
    
    mw=Chem.RWMol(m)
    if substruct==None:
        
        #check if it has any of the substructures in the atom list
        no_match=True 
        print(list(md.substruct_dict.keys()))
        for s in list(md.substruct_dict.keys()):
            s_mol=md.substruct_dict[s]
            print(mw.GetSubstructMatches(s_mol))
            for match in mw.GetSubstructMatches(s_mol):
                if match[0] in at_list:
                    no_match=False
                    break
                if not no_match:
                    break
        
        #if there is at least one substructure match in the molecule within the atom list:
        if no_match==False:
            
            #assume the current substructure cannot be found within the atom list
            no_match_in_list=True
            
            #while the current substructure is not in the atom list, pick a random substructure
            while no_match_in_list:
                substruct=choice(list(md.substruct_dict.keys()))
                substr_mol=md.substruct_dict[substruct]
                for match in mw.GetSubstructMatches(substr_mol):
                    if match[0] in at_list:
                        no_match_in_list=False
        else:
            print('no match') #if it doesn't have any of the substructures
            return None

    else: #if the substructure is specified
        substr_mol=md.substruct_dict[substruct]
    
    #find a match in the list of atoms
    no_match_in_list=True
    a=mw.GetSubstructMatches(substr_mol)
    
    for b in a:
        if b[0] in at_list:
            no_match_in_list=False
            break
    
    if not no_match_in_list:
        at_idx=choice(a)[0]
        while at_idx not in at_list:
            at_idx=choice(a)[0]
    else:
        return None #if the structure cannot be found within the atom list
    
    #pick a random modification function for the given substructure
    fct=choice(md.modif_dict[substr_mol])
    new_m=fct(mw,at_idx)
    
    mw_mfp=AllChem.GetMorganFingerprint(mw,4)
    new_m_mfp=AllChem.GetMorganFingerprint(new_m,4)
    
    #while the similarity between the new structure and the original one is lower than the similarity threshold,
    #get a new modified structure
    
    j=0
    while DataStructs.DiceSimilarity(mw_mfp,new_m_mfp)<similarity and j<50:
        j=j+1
        fct=choice(md.modif_dict[substr_mol])
        new_m=fct(mw,at_idx)
        
        mw_mfp=AllChem.GetMorganFingerprint(mw,4)
        new_m_mfp=AllChem.GetMorganFingerprint(new_m,4)
    
    if DataStructs.DiceSimilarity(mw_mfp,new_m_mfp)>=similarity:
        return new_m
    elif j>=50:
        return None

        
def weighted(m,at_list):
    """
    arguments: mol object, list of atom indices
    
    returns: new mol object or None (if new mol object not possible)
        
    this function performs a replace_substituent or CHtoN modification on the molecule (the change is
    picked randomly from a weighted list of changes)
    """
    
    mw=Chem.RWMol(m)
    no_match=True   #check if there is any modification possible
    for s in ['cH','c']:
        s_mol=md.substruct_dict[s]
        for match in mw.GetSubstructMatches(s_mol):
            if match[0] in at_list:
                no_match=False
                break
            if not no_match:
                break
    if no_match:
        print('no match')
        return None
    else:
        
        #find a match in the list of atoms
        
        no_match_in_list=True
        new_m=None
        
        while no_match_in_list:
            change=choice(weighted_list)
            if change=='-n':
                
                substr_mol=md.substruct_dict['cH']
        
                a=mw.GetSubstructMatches(substr_mol)
                
                for b in a:
                    if b[0] in at_list:
                        no_match_in_list=False
                        break
                
                if not no_match_in_list:
                    print('chose '+change)
                    at_idx=choice(a)[0]
                    while at_idx not in at_list:
                        at_idx=choice(a)[0]
        
                    new_m=md.CHToN(mw,at_idx)
                
            else:
                substr_mol=md.substruct_dict['c']
        
                a=mw.GetSubstructMatches(substr_mol)
                
                for b in a:
                    if b[0] in at_list:
                        no_match_in_list=False
                        break
                
                if not no_match_in_list:
                    print('chose '+change)
                    at_idx=choice(a)[0]
                    while at_idx not in at_list:
                        at_idx=choice(a)[0]
        
                    new_m=md.replace_substituent(mw,at_idx,subst=change)
                    
        return new_m
          
        

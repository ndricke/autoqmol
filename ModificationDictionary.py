#!/usr/bin/env python

"""
arguments: -
    
calls:
    
this script contains functions for modifying mol objects and a dictionary of substructures to be found in the 
mol object and modification functions associated with each substructure

!!! The list of atomic numbers to be used when replacing H with another atom has to be changed manually !!!

"""
from rdkit import Chem
from rdkit.Chem import AllChem
from random import choice
from collections import OrderedDict

#a list of atomic numbers to be used when replacing H with another atom
#Z_list=[5,6,7,8,9,14,15,16,17,35,53]
Z_list=[6,7,8,9,17,35,53]
Mod_list=['CN','OCH3','CH3','NH2','OH','F','Cl','Br','I']

#a dictionary of functional groups to be used when modifying structures
fg=OrderedDict()
fName='/home/nricke/ulysses/iandriuc/script_withoutscaffold/FunctionalGroupsForModifying.txt'
with open(fName) as f:
    for line in f:
        fg[line.split()[0]]=Chem.MolFromSmarts(line.split()[1])


def HToOtherElement(m,cn_idx,Z=None):
    """
    arguments: mol object, connection atom index, Z (optional)
    
    returns: new mol object
    
    this function replaces a H atom bound to an atom with index cn_idx with another atom with atomic number Z
    """
    print('chose HToOtherElement')
    if Z==None:
        Z=choice(Z_list)  #pick random element
    mw=Chem.RWMol(m)
    for at in mw.GetAtomWithIdx(cn_idx).GetNeighbors():
        if at.GetSymbol()=='H':
            H_idx=at.GetIdx()
            break
        
    mw.ReplaceAtom(H_idx,Chem.Atom(Z))
    Chem.SanitizeMol(mw)
    
    
    #add Hs to satisfy the valence of the new atom
    while len(mw.GetAtomWithIdx(H_idx).GetNeighbors())<mw.GetAtomWithIdx(H_idx).GetTotalValence():
        idx=mw.AddAtom(Chem.Atom(1))
        mw.AddBond(H_idx,idx,Chem.BondType.SINGLE)

    Chem.SanitizeMol(mw)  
    AllChem.EmbedMolecule(mw)
    AllChem.MMFFOptimizeMolecule(mw)
      
    return mw
    
def HToPh(m,cn_idx):
    """
    arguments: mol object, connection atom index 
    
    returns: new mol object
    
    this function replaces a H atom bound to an atom with index cn_idx with a phenyl ring
    """
    print('chose HtoPh')

    mw=Chem.RWMol(m)
    
    #remove the H atom
    for at in mw.GetAtomWithIdx(cn_idx).GetNeighbors():
        if at.GetSymbol()=='H':
            mw.RemoveAtom(at.GetIdx())

    #build the Ph ring        
    c1=mw.AddAtom(Chem.Atom(6))
    mw.AddBond(c1,cn_idx,Chem.BondType.SINGLE)
    
    c2=mw.AddAtom(Chem.Atom(6))
    mw.AddBond(c2,c1,Chem.BondType.DOUBLE)
    h2=mw.AddAtom(Chem.Atom(1))
    mw.AddBond(h2,c2,Chem.BondType.SINGLE)
    
    c3=mw.AddAtom(Chem.Atom(6))
    mw.AddBond(c2,c3,Chem.BondType.SINGLE)
    h3=mw.AddAtom(Chem.Atom(1))
    mw.AddBond(h3,c3,Chem.BondType.SINGLE)
    
    c4=mw.AddAtom(Chem.Atom(6))
    mw.AddBond(c4,c3,Chem.BondType.DOUBLE)
    h4=mw.AddAtom(Chem.Atom(1))
    mw.AddBond(h4,c4,Chem.BondType.SINGLE)    

    c5=mw.AddAtom(Chem.Atom(6))
    mw.AddBond(c5,c4,Chem.BondType.SINGLE)
    h5=mw.AddAtom(Chem.Atom(1))
    mw.AddBond(h5,c5,Chem.BondType.SINGLE)

    c6=mw.AddAtom(Chem.Atom(6))
    mw.AddBond(c6,c5,Chem.BondType.DOUBLE)
    h6=mw.AddAtom(Chem.Atom(1))
    mw.AddBond(h6,c6,Chem.BondType.SINGLE) 
    
    mw.AddBond(c6,c1,Chem.BondType.SINGLE)

    frozen_dict={}

    for j in range(len(m.GetAtoms())):
        if j<c1 and j!=cn_idx:
            frozen_dict[j]=m.GetConformer(0).GetAtomPosition(j)
            
    print(Chem.MolToSmiles(mw))

    Chem.SanitizeMol(mw)  
    AllChem.EmbedMolecule(mw,coordMap=frozen_dict)

    return mw    

def CHToN(m,cn_idx,subst=None):
    
    """
    arguments: mol object, connection atom index, substituent (optional and not necessary for this function,
    but when called from Modifications2 it requires subst as an argument)
    
    returns: new mol object
    
    this function replaces an aromatic atom with index cn_idx with a N atom bound to a substituent
    """
    
    print('chose CHtoN')
    mw=Chem.RWMol(m)
    Chem.Kekulize(mw,clearAromaticFlags=True)
    #the ring the atom is a part of 
    for ring in mw.GetRingInfo().AtomRings():
        if cn_idx in ring:
            target_ring=ring   #ring containing the connection atom
   
    #remove the ngbs that are not part of the ring
    for ngb in mw.GetAtomWithIdx(cn_idx).GetNeighbors():
        if ngb.GetIdx() not in target_ring:
            mw.RemoveAtom(ngb.GetIdx())
    
    #replace C with N
    mw.ReplaceAtom(cn_idx,Chem.Atom(7))
#    mw.GetAtomWithIdx(cn_idx).SetFormalCharge(1)
#    h=mw.AddAtom(Chem.Atom(1))
#    mw.AddBond(h,cn_idx,Chem.BondType.SINGLE)
#    Chem.SanitizeMol(mw)
    
    return mw

def CHToNplus(m,cn_idx,subst=None):
    
    """
    arguments: mol object, connection atom index, substituent (optional)
    
    returns: new mol object
    
    this function replaces an aromatic atom with index cn_idx with a N atom bound to a substituent
    """
    
    print('chose CHtoNplus')
    mw=Chem.RWMol(m)
    Chem.Kekulize(mw,clearAromaticFlags=True)
    #the ring the atom is a part of 
    for ring in mw.GetRingInfo().AtomRings():
        if cn_idx in ring:
            target_ring=ring
   
    #remove the ngbs that are not part of the ring
    for ngb in mw.GetAtomWithIdx(cn_idx).GetNeighbors():
        if ngb.GetIdx() not in target_ring:
            mw.RemoveAtom(ngb.GetIdx())
    
    #replace C with N
    mw.ReplaceAtom(cn_idx,Chem.Atom(7))
    mw.GetAtomWithIdx(cn_idx).SetFormalCharge(1)
    h=mw.AddAtom(Chem.Atom(1))
    mw.AddBond(h,cn_idx,Chem.BondType.SINGLE)
#    Chem.SanitizeMol(mw)
    
    if subst=='H':
        Chem.SanitizeMol(mw)  
        AllChem.EmbedMolecule(mw)

        return mw
    elif subst==None:
        fct=choice(modif_dict[substruct_dict['n+H']])
        if fct==HToOtherElement:
            mw=fct(mw,cn_idx,6)
        else:
            mw=fct(mw,cn_idx)

    
          
        return mw
    else:
        print('This substituent is not supported yet')
    
def AddSubstToPh(m,cn_idx,Z=None,pos=None):

    """
    arguments: mol object, connection atom index, Z (optional), pos (optional - o/m/p)
    
    returns: new mol object
    
    this function replaces an H atom at position pos on a phenyl ring bound to an atom with index cn_idx with
    another atom with atomic number Z
    """
    if Z==None:
        Z=choice(Z_list)
    if pos==None:
        pos=choice(['o','m','p'])
        
    print('chose AddSubstToPh at position '+pos)
    
    mw=Chem.RWMol(m)
    
    #the ring the atom is a part of 
    for ring in mw.GetRingInfo().AtomRings():
        if cn_idx in ring:
            target_ring=ring

    
    for at in mw.GetAtomWithIdx(cn_idx).GetNeighbors():
        if at.GetIdx() not in target_ring:
            Ph_C_idx=at.GetIdx()    #The C atom that is part of the Ph ring and is bound to the connection atom
            break
        
    for ngb1 in mw.GetAtomWithIdx(Ph_C_idx).GetNeighbors():
        if ngb1.GetSymbol()=='C' and ngb1.GetIdx()!=cn_idx:
            o_c=ngb1.GetIdx()   #the carbon ortho to Ph_C_idx
    for ngb2 in mw.GetAtomWithIdx(o_c).GetNeighbors():
        if ngb2.GetSymbol()=='H':
            o_h=ngb2.GetIdx()   #the hydrogen ortho to Ph_C_idx
        elif ngb2.GetIdx()!=Ph_C_idx:
            m_c=ngb2.GetIdx()   #the carbon meta to Ph_C_idx
    for ngb3 in mw.GetAtomWithIdx(m_c).GetNeighbors():
        if ngb3.GetSymbol()=='H':
            m_h=ngb3.GetIdx()   #the hydrogen meta to Ph_C_idx
        elif ngb3.GetIdx()!=o_c:
            p_c=ngb3.GetIdx()   #the carbon para to Ph_C_idx
    for ngb4 in mw.GetAtomWithIdx(p_c).GetNeighbors():
        if ngb4.GetSymbol()=='H':
            p_h=ngb4.GetIdx()   #the hydrogen para to Ph_C_idx
                    
    if pos=='o':
        H_idx=o_h
    elif pos=='m':
        H_idx=m_h
    elif pos=='p':
        H_idx=p_h
        
    mw.ReplaceAtom(H_idx,Chem.Atom(Z))
    Chem.SanitizeMol(mw)
    
    
    #add Hs to satisfy the valence of the new atom
    while len(mw.GetAtomWithIdx(H_idx).GetNeighbors())<mw.GetAtomWithIdx(H_idx).GetTotalValence():
        idx=mw.AddAtom(Chem.Atom(1))
        mw.AddBond(H_idx,idx,Chem.BondType.SINGLE)

    Chem.SanitizeMol(mw)  
    AllChem.EmbedMolecule(mw)
    AllChem.MMFFOptimizeMolecule(mw)
      
    return mw

def replace_substituent(m,cn_idx,subst=None):
    """
    arguments: mol object, connection atom index, substituent (optional)
        
    returns: new mol object
        
    this function replaces a substituent on an aromatic carbon with the given substituent or a random one
    """
    print('chose replace_substituent')
    mw=Chem.RWMol(m)
#    print('cn idx: '+str(cn_idx))
    for ring in mw.GetRingInfo().AtomRings():
        if cn_idx in ring:
            target_ring=ring #the ring the connection atom is part of
            break
        
    for ngb in mw.GetAtomWithIdx(cn_idx).GetNeighbors():
        if ngb.GetIdx() not in target_ring:
            subst_idx=ngb.GetIdx()

    reached_the_end=False
#    print('subst idx: '+str(subst_idx))
    current_atom_list=[]
    all_atoms=[subst_idx]
    for at in mw.GetAtomWithIdx(subst_idx).GetNeighbors():
        if at.GetIdx()!=cn_idx:
            current_atom_list.append(at.GetIdx())
#    print('current atom list: ')
#    print([str(x) for x in current_atom_list])
    all_atoms.extend(current_atom_list)     #a list of all the atoms making up the substituent to be replaced
#    print('all atoms: ')
#    print([str(x) for x in all_atoms])

    while not reached_the_end:
        new_atom_list=[]
        for at_idx in current_atom_list:
            
            for n in mw.GetAtomWithIdx(at_idx).GetNeighbors():
                if n.GetIdx() not in all_atoms and n.GetIdx()!=subst_idx and n.GetIdx() not in new_atom_list:
                    new_atom_list.append(n.GetIdx())
                    
#        print('new atom list: ')
#        print([str(x) for x in new_atom_list])
        if len(new_atom_list)==0:
            reached_the_end=True
        else:
            
            current_atom_list=new_atom_list
#            print('current atom list: ')
#            print([str(x) for x in current_atom_list])
            all_atoms.extend(current_atom_list)
#            print('all atoms: ')
#            print([str(x) for x in all_atoms])
    
    new_cn_idx=cn_idx
    all_atoms.sort(reverse=True) #remove them in descending order so that the indices of the atoms that have
                                 #yet to be removed does not change
    for at_idx in all_atoms:
        mw.RemoveAtom(at_idx)

        if at_idx<cn_idx:
            new_cn_idx=new_cn_idx-1
        
    if subst==None:
        k=choice(fg.keys())
        print('chose '+k)
        #the idx of the connection atom from the incoming FG
        i=len(mw.GetAtoms())
        
        a=Chem.RWMol(fg[k])
        a.RemoveAtom(0)
        new_m=Chem.CombineMols(mw,a)    #add new FG
        new_m=Chem.RWMol(new_m)
        new_m.AddBond(new_cn_idx,i,Chem.BondType.SINGLE)
#        print(Chem.MolToSmiles(new_m))
        Chem.SanitizeMol(new_m)  
        AllChem.EmbedMolecule(new_m)
        AllChem.MMFFOptimizeMolecule(new_m)

        return new_m
    elif '-'+subst or subst in fg.keys():
        if subst[0]=='-':
            k=subst
        else:
            k='-'+subst
        
        #the idx of the connection atom from the incoming FG
        i=len(mw.GetAtoms())
        
        a=Chem.RWMol(fg[k])
        a.RemoveAtom(0)
        new_m=Chem.CombineMols(mw,a)
        new_m=Chem.RWMol(new_m)
        new_m.AddBond(new_cn_idx,i,Chem.BondType.SINGLE)
        Chem.SanitizeMol(new_m)  
        AllChem.EmbedMolecule(new_m)
        AllChem.MMFFOptimizeMolecule(new_m)

        return new_m

    else:
        print('This substituent is not supported yet')


# a dictionary which maps mol objects onto strings representing substructures
substruct_dict={}

# a dictionary in which each substructure mol object from substruct_dict is assigned a list of possible modification
# functions
modif_dict={}

#substruct_dict['cH']=Chem.MolFromSmarts('[c;H1,C;H1]')
#substruct_dict['cH']=Chem.MolFromSmarts('[c,C;H1]')
substruct_dict['cH']=Chem.MolFromSmarts('[c;H1]')
modif_dict[substruct_dict['cH']]=[]
#modif_dict[substruct_dict['cH']].append(HToOtherElement)
modif_dict[substruct_dict['cH']].append(replace_substituent)
#modif_dict[substruct_dict['cH']].append(CHToN)

#substruct_dict['cH']=Chem.MolFromSmarts('[cX3H]')
#modif_dict[substruct_dict['cH']]=[]
#modif_dict[substruct_dict['cH']].append(HToOtherElement)
#modif_dict[substruct_dict['cH']].append(CHToN)
#modif_dict[substruct_dict['cH']].append(CHToNplus)

#substruct_dict['n+H']=Chem.MolFromSmarts('[n+H]')
#modif_dict[substruct_dict['n+H']]=[]
#modif_dict[substruct_dict['n+H']].append(HToOtherElement)
#modif_dict[substruct_dict['n+H']].append(HToPh)
#
#substruct_dict['n+Ph']=Chem.MolFromSmarts('[n+]-[cX3]1:[cX3H]:[cX3H]:[cX3H]:[cX3H]:[cX3H]1')
#modif_dict[substruct_dict['n+Ph']]=[]
#modif_dict[substruct_dict['n+Ph']].append(AddSubstToPh)
#
#substruct_dict['c']=Chem.MolFromSmarts('[c;R1;!$(c=O);!$(c=S)]')
#modif_dict[substruct_dict['c']]=[]
#modif_dict[substruct_dict['c']].append(replace_substituent)

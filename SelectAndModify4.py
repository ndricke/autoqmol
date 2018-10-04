#!/usr/bin/env python
"""
arguments: parent directory, destination directory, number of modifications to be made for each catalyst, which 
method to consider when selecting which catalysts bind O2 (b3lyp/sp/random, optional, random by default), number of
random catalysts to be selected (if the selection method is 'random'), which method to use when modifying the
catalysts (random/weighted/similarity, optional, random by default), similarity index (if the modification method
is 'similarity', optional, between 0 and 1, 0 by default)
    
calls: Modifications2, add_o2_frozen, multi_jobs
    
this script selects the catalysts that were predicted to bind O2 by one method and modifies them a given number
of times (i.e. creates a number of new catalysts equal to 'number of modifications')
it creates directories, mol and xyz files for each catalyst and O2 bound species and runs optimizations (b3lyp,
3-21g, no solvation)

if selection is not random, the modification is made on target (i.e. around the binding site)

if selection is random, the modifcation is made anywhere outside the scaffold

!!! The cutoff index for the scaffold has to be changed manually !!!
"""
import os
import pandas as pd
from rdkit import Chem
from Modifications2 import random_on_target
from Modifications2 import random_from_list
from Modifications2 import weighted
from shutil import copy
import argparse
from random import choice

parser = argparse.ArgumentParser()
parser.add_argument("parent_dir", type=str, help="the parent directory from which to select catalysts to be modified")
parser.add_argument("mod_parent_dir", type=str, help="the name of the destination directory for the modified catalysts")
parser.add_argument("modif_no", type=int, help="number of modifications to be applied to each catalyst")
parser.add_argument("selection",type=str, choices=['b3lyp','sp','random'],nargs='?',default='random',help="method for selecting catalysts to be modified (random by default)")
parser.add_argument("-s",type=int,default=None,help="number of random selections to be made (if the selection method is random)")
parser.add_argument("modification",type=str,choices=['weighted','similarity','random'],nargs='?',default='random',help="method for selecting what modifications to make (random by default)")
parser.add_argument("-i",type=float,default=None,help="similarity index (if the modification method is 'similarity')")
args = parser.parse_args()

# !!! CHANGE MANUALLY !!!
cutoff=7

if args.selection=="random" and args.s==None:
    print('Please specify the number of catalysts to be randomly selected')
elif args.modification=='similarity' and args.i==None:
    print('Please specify the similarity index')
elif args.i!=None and (args.i<=0 or args.i>=1):
    print('The similarity index has to be between 0 and 1')
else: 
    print(args)

    parent_dir=args.parent_dir
    
    modif_no=args.modif_no
    
    method=args.selection
    
    modif=args.modification
    
    if args.mod_parent_dir[-1]=='/':
        mod_parent_dir=args.mod_parent_dir[:-1]
    else:
        mod_parent_dir=args.mod_parent_dir
    
    try:
        os.mkdir(os.getcwd()+'/'+mod_parent_dir)
    except:
        pass
    
    #copy O2 files to the new directory
    for f in os.listdir(parent_dir):
        if 'o2' in f and '_out' not in f:
            copy(parent_dir+'/'+f,mod_parent_dir+'/')
            
    if method!='random':
        

        sum_file_path=parent_dir+"/summary.csv"
        
        #read dataframe from summary file
        df=pd.read_csv(sum_file_path,sep='\t',index_col=None)
        
        for col_name in list(df.columns.values):
            if "Y/N" in col_name and method in col_name:
                tgt_col_name=col_name
                break
            
        for i in range(len(df)):
            
            if df.loc[i][tgt_col_name].strip()=="y":    #if it is predicted to bind O2
                
                c_name=df.loc[i][0].strip()     #candidate name
                in_site=str(df.loc[i][1]).strip()            #initial binding site
                candidate_dir=parent_dir+"/"+c_name
                bare_dir=candidate_dir+'/bare'
                o2_dir=candidate_dir+'/o2'
                
                c_mol_file=bare_dir+'/'+c_name.replace('c_','cat')+'.mol'
                m=Chem.MolFromMolFile(c_mol_file, removeHs=False)
                
                o2_out_mol_file=o2_dir+'/'+c_name.replace('c_','cat')+'_'+in_site+'_out.mol'
                m_o2=Chem.MolFromMolFile(o2_out_mol_file, removeHs=False)
                
                mod_list=[Chem.MolToSmiles(m)]      #list of canonical smiles of the modified structures   
                at_list=list(range(cutoff+1,len(m.GetAtoms())))
                if modif=='random':
                    n=random_on_target(m,m_o2,cutoff)
                elif modif=='similarity':
                    n=random_from_list(m,at_list,similarity=args.i)
                elif modif=='weighted':
                    n=weighted(m,at_list)


                
                for j in range(modif_no):
                    
                    print('working on modif '+str(j+1))
                    
                    if n!=None:
                    
                        while Chem.MolToSmiles(n) in mod_list: #if already in list, try a new one
                            if modif=='random':
                                n=random_on_target(m,m_o2,cutoff)
                            elif modif=='similarity':
                                n=random_from_list(m,at_list,similarity=args.i)
                            elif modif=='weighted':
                                n=weighted(m,at_list)
                        mod_list.append(Chem.MolToSmiles(n))
                        
                        mod_c_dir=mod_parent_dir+'/'+c_name+'_m'+(2-len(str(j+1)))*'0'+str(j+1)
                        mod_bare_dir=mod_c_dir+'/bare'
                        mod_o2_dir=mod_c_dir+'/o2'
                        
                        xyz_file=mod_bare_dir+'/'+c_name.replace('c_','cat')+'_m'+(2-len(str(j+1)))*'0'+str(j+1)+'.xyz'
                        mol_file=mod_bare_dir+'/'+c_name.replace('c_','cat')+'_m'+(2-len(str(j+1)))*'0'+str(j+1)+'.mol'
                        
                        try:
                            os.mkdir(mod_c_dir)
                            os.mkdir(mod_bare_dir)
                            os.mkdir(mod_o2_dir)
                        except:
                            pass
                        
                        #write mol file for bare modified catalyst
                        with open(mol_file,'w') as g:
                            g.write(Chem.MolToMolBlock(n))
                            
                        #mol to xyz
                        os.system("babel -imol %s -oxyz %s" % (mol_file,xyz_file))
        #                print("babel -imol %s -oxyz %s" % (mol_file,xyz_file))
                        
                        #run add_o2_frozen which creates .xyz and .mol files for the bare and o2 bound structures
                        os.system('python add_o2_frozen.py '+mod_c_dir+' '+str(cutoff))
        #                print('python add_o2_frozen.py '+mod_c_dir+' '+str(cutoff))
                        
                        #run optimisations for the bare structure
                        os.system('python multi_jobs.py '+mod_bare_dir+' b3lyp 3-21g opt - n')
        #                print('python multi_jobs.py '+mod_bare_dir+' b3lyp 3-21g opt - n')
                        
                        #run optimisations for the o2 bound structures
                        os.system('python multi_jobs.py '+mod_o2_dir+' b3lyp 3-21g opt - n')
        #                print('python multi_jobs.py '+mod_o2_dir+' b3lyp 3-21g opt - n')
        #                print('\n\n')
                        
                    else:
                        print('cannot modify '+c_name+' anymore')
                    
    elif method=='random':
        sum_file_path=parent_dir+"/summary.csv"
        
        #read dataframe from summary file
        df=pd.read_csv(sum_file_path,sep='\t',index_col=None)
        
        cat_list=[]
        for i in range(len(df)):
            
            cat_list.append(df.loc[i][0].strip())
            
        random_list=[]
        rand_cat=choice(cat_list)
        for i in range(args.s):
            while rand_cat in random_list:
                rand_cat=choice(cat_list)
            random_list.append(rand_cat)
        
        
        for c_name in random_list:
            print('\n'+c_name+'\n\n')
            candidate_dir=parent_dir+"/"+c_name
            bare_dir=candidate_dir+'/bare'
            o2_dir=candidate_dir+'/o2'
            
            c_mol_file=bare_dir+'/'+c_name.replace('c_','cat')+'.mol'
            m=Chem.MolFromMolFile(c_mol_file, removeHs=False)

            mod_list=[Chem.MolToSmiles(m)]     

            at_list=list(range(cutoff+1,len(m.GetAtoms())))
    
            if modif=='random':
                n=random_on_target(m,m_o2,cutoff)
            elif modif=='similarity':
                n=random_from_list(m,at_list,similarity=args.i)
            elif modif=='weighted':
                n=weighted(m,at_list)
            
            for j in range(modif_no):
                
                print('working on modif '+str(j+1))
                
                if n!=None:
                
                    while Chem.MolToSmiles(n) in mod_list:
                        if modif=='random':
                            n=random_on_target(m,m_o2,cutoff)
                        elif modif=='similarity':
                            n=random_from_list(m,at_list,similarity=args.i)
                        elif modif=='weighted':
                            n=weighted(m,at_list)
                        
                    mod_list.append(Chem.MolToSmiles(n))
                    
                    mod_c_dir=mod_parent_dir+'/'+c_name+'_m'+(2-len(str(j+1)))*'0'+str(j+1)
                    mod_bare_dir=mod_c_dir+'/bare'
                    mod_o2_dir=mod_c_dir+'/o2'
                    
                    xyz_file=mod_bare_dir+'/'+c_name.replace('c_','cat')+'_m'+(2-len(str(j+1)))*'0'+str(j+1)+'.xyz'
                    mol_file=mod_bare_dir+'/'+c_name.replace('c_','cat')+'_m'+(2-len(str(j+1)))*'0'+str(j+1)+'.mol'
                    
                    try:
                        os.mkdir(mod_c_dir)
                        os.mkdir(mod_bare_dir)
                        os.mkdir(mod_o2_dir)
                    except:
                        pass
                    
                    #write mol file for bare modified catalyst
                    with open(mol_file,'w') as g:
                        g.write(Chem.MolToMolBlock(n))
                        
                    #mol to xyz
                    os.system("babel -imol %s -oxyz %s" % (mol_file,xyz_file))
    #                print("babel -imol %s -oxyz %s" % (mol_file,xyz_file))
                    
                    #run add_o2_frozen which creates .xyz and .mol files for the bare and o2 bound structures
                    os.system('python add_o2_frozen.py '+mod_c_dir+' '+str(cutoff))
    #                print('python add_o2_frozen.py '+mod_c_dir+' '+str(cutoff))
                    
                    #run optimisations for the bare structure
                    os.system('python multi_jobs.py '+mod_bare_dir+' b3lyp 3-21g opt - n')
    #                print('python multi_jobs.py '+mod_bare_dir+' b3lyp 3-21g opt - n')
                    
                    #run optimisations for the o2 bound structures
                    os.system('python multi_jobs.py '+mod_o2_dir+' b3lyp 3-21g opt - n')
    #                print('python multi_jobs.py '+mod_o2_dir+' b3lyp 3-21g opt - n')
    #                print('\n\n')
                    
                else:
                    print('cannot modify '+c_name+' anymore')


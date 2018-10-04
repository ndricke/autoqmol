#!/usr/bin/env python2.7
import pyQChem as qc
import numpy as np
import sys
import os
import copy


charge_tran = {'a2':'-2','a1':'-1','a0':'0','c1':'1','c2':'2','c3':'3','c4':'4'}
revd=dict([reversed(i) for i in charge_tran.items()])
charge_tran.update(revd)

class QcIn(object):

    def __init__(self, infile,charge=0,mult=1,jobtype='opt',basis='6-31+g*',method='tpssh', \
                 geom_read=False, nametrunc=True, bond_cons = None, infodump=False):

        self.jobdict = {'sp':self.jobSp, 'opt':self.jobOpt, 'fq':self.jobFreq, 'sopt':self.jobSopt, \
                        'cube':self.jobCube, 'constr':self.jobConstrainedOpt, 'fixbonds':self.jobFixedBondOpt, \
                        'cdft':self.jobCdft, \
                       }

#        self.bond_cons = [[2,5,con_N_x],[3,4,con_N_x],[2,3,con_N_y],[4,5,con_N_y]]
        self.bond_cons = bond_cons
        self.info_dump = infodump

#        self.free_atom = ['Fe','O']
        self.free_atom = ['O']
        self.job_arr_list = []

        self.charge = charge
        self.mult = mult
        self.jobtype = jobtype
        self.pcm_arr = None; self.sol_arr = None; self.plot_arr = None

        #if nametrunc == True: self.name = (infile.split('.')[0]).split('_')[0]
        #else: self.name = infile.split('.')[0]

        splinfile = infile.split('/')
        if len(splinfile) == 1:
            self.path = ''
        else:
            self.path = '/'.join(splinfile[:-1])+'/'

        if nametrunc == True: self.name = (splinfile[-1].split('.')[0]).split('_')[0]
        else: self.name = (splinfile[-1].split('.'))[0]


        print(splinfile)
        print(splinfile[-1])
        print(self.name)

        print(self.charge)
        print(self.mult)
        self.out = qc.read(infile)
        self.rem = qc.rem_array()
        self.rem.basis(basis)
        self.rem.method(method)
        self.rem.thresh("14")
        self.rem.add('mem_total','4096')
        self.rem.add('mem_static','256')
        self.rem.add('max_scf_cycles','500')
#        self.rem.add('scf_convergence','7')
#        self.rem.add('scf_algorithm','diis_gdm')

        ##For these metal-centered systems
        self.rem.add('unrestricted','true')
#        if self.mult == '1': ## I didn't find that this worked at all for any of the metal-macrocycle systems
#            self.rem.add("scf_guess_mix", "1")
#            self.rem.add("scf_guess", "gwh")

        if geom_read:
            self.mol = qc.mol_array()
            self.mol.geometry("read")
            self.rem.add('scf_guess', 'read')
        else:
            if infile.split('.')[-1] == 'out':
#                self.mol=qc.mol_array(self.out.opt.geometries[-1])
                self.mol = qc.mol_array(self.out.general.final_geometry)
#                self.mol = qc.mol_array(self.out.general.initial_geometry)
            else:
                xyz = qc.cartesian(atom_list=self.out.list_of_atoms)
                self.mol = qc.mol_array(xyz)
        self.mol.charge(charge)
        self.mol.multiplicity(mult)

    def addSolvent(self):
        self.rem.solvent_method('PCM')
        self.pcm_arr = qc._unsupported_array("pcm")
        self.pcm_arr.add_line("theory cpcm")
        self.pcm_arr.add_line("method swig")
        self.pcm_arr.add_line("solver inversion")
        self.pcm_arr.add_line("heavypoints 194")
        self.pcm_arr.add_line("hpoints 194")
        self.pcm_arr.add_line("radii bondi")
        self.pcm_arr.add_line("vdwscale 1.2")

        self.sol_arr = qc._unsupported_array("solvent")
        self.sol_arr.add_line("dielectric 78.39")

        self.vdw_arr = qc._unsupported_array("van_der_waals")
        self.vdw_arr.add_line("1")

        """Source: http://periodictable.com/Properties/A/VanDerWaalsRadius.v.html"""
        #self.vdw_arr.add_line("24  1.97")
        #self.vdw_arr.add_line("25  1.96")
        #self.vdw_arr.add_line("26  1.96")
        #self.vdw_arr.add_line("27  1.95")
        #self.vdw_arr.add_line("28  1.63")
        #self.vdw_arr.add_line("29  1.40")
        #self.vdw_arr.add_line("30  1.39")

        """Source: https://physlab.lums.edu.pk/images/f/f6/Franck_ref2.pdf0""Source: http://periodictable.com/Properties/A/VanDerWaalsRadius.v.html"""
        self.vdw_arr.add_line("24  1.97")
        self.vdw_arr.add_line("25  1.96")
        self.vdw_arr.add_line("26  1.96")
        self.vdw_arr.add_line("27  1.95")
        self.vdw_arr.add_line("28  1.94")
        self.vdw_arr.add_line("29  2.00")
        self.vdw_arr.add_line("30  2.02")

        self.job_arr_list.append(self.pcm_arr)
        self.job_arr_list.append(self.sol_arr)
        self.job_arr_list.append(self.vdw_arr)

    @staticmethod
    def splitCoord(coord_list):
        coords = np.zeros((len(coord_list),3))
        atom_list = []
        for i, item in enumerate(coord_list):
            atom_list.append(item[0])
            coords[i,:] = [float(c) for c in item[1:]]
        return atom_list, coords

    def O2coord(self):
        atom_list, atom_coord = self.splitCoord(self.out.list_of_atoms)
        O_indices = [i for i, x in enumerate(atom_list) if x[0] == 'O']
        if len(O_indices) >= 2:
            for index1 in range(len(O_indices)-1):
                for index2 in range(index1+1, len(O_indices)):
                    if abs(np.linalg.norm(atom_coord[O_indices[index1],:] - atom_coord[O_indices[index2],:])) < 1.6:
                        return O_indices[index1], O_indices[index2]
        else:
            print("Couldn't find the O2.")
            return None

#CDFT calculation
    def jobCdft(self):
        self.rem.add('CDFT', 'True')
        self.rem.add('SCF_PRINT', 'True')
        self.cdft_arr = qc._unsupported_array("cdft")
        self.cdft_arr.add_line("1")
        O2_indices = self.O2coord()
        #print("O2 indices:", O2_indices)
        self.cdft_arr.add_line("1  %s  %s  s" % (O2_indices[0]+1, O2_indices[1]+1))
        # error if O2coord() returns None
        self.job_arr_list.append(self.cdft_arr)

        atoms = copy.copy(self.out.list_of_atoms)
        first_two = [atoms[0][0], atoms[1][0]]
        last_two = [atoms[-2][0], atoms[-1][0]]
        #if first_two != ['O','O'] and last_two == ['O','O']:
        #    atoms.insert(0, atoms.pop()); atoms.insert(0, atoms.pop())
        #xyz = qc.cartesian(atom_list=atoms)
        #self.mol = qc.mol_array(xyz)
        #self.mol.charge(self.charge)
        #self.mol.multiplicity(self.mult)

#Fix all but selected atom numbers
    def jobConstrainedOpt(self):
        self.rem.jobtype("opt")
        self.con_arr = qc._unsupported_array("opt")
        self.con_arr.add_line("FIXED")
        for i,atom in enumerate(self.out.list_of_atoms):
            if atom[0] in self.free_atom:
                continue
            else:
                self.con_arr.add_line(str(i+1) + " XYZ")
        self.con_arr.add_line("ENDFIXED")
        self.job_arr_list.append(self.con_arr)

#Fix bond lengths
    def jobFixedBondOpt(self):
        self.rem.jobtype("opt")
        self.con_arr = qc._unsupported_array("opt")
        self.con_arr.add_line("CONSTRAINT")
        for constr in self.bond_cons: #need to implement this to read these from file if we used this for other systems
            self.con_arr.add_line("stre %i %i %f" % (constr[0], constr[1], constr[2]))
        self.con_arr.add_line("ENDCONSTRAINT")
        self.job_arr_list.append(self.con_arr)

    def jobCube(self):
        self.jobSp()
        self.rem.add("make_cube_files", "true")
        self.plot_arr = qc._unsupported_array("plots")
        self.plot_arr.add_line("Grid information comment")
        self.plot_arr.add_line("200 -10.0 10.0")
        self.plot_arr.add_line("200 -10.0 10.0")
        self.plot_arr.add_line("200 -10.0 10.0")
        self.plot_arr.add_line("1 0 0 0")
        self.plot_arr.add_line("1")
        self.job_arr_list.append(self.plot_arr)

    #Single point job
    def jobSp(self):
        self.rem.jobtype("sp")
        self.rem.add("print_orbitals","true")
        self.rem.add("chelpg", "true")
        #self.rem.add("lowdin_population", "true")

        if self.info_dump == True:
            self.rem.add("gui=2", "")
            self.rem.add("molden_format", "true")

        self.addSolvent()

    #Geometry Optimization Job
    def jobOpt(self):
        self.rem.jobtype("opt")
#        self.rem.add("geom_opt_tol_gradient", "5")
#        self.rem.add("geom_opt_tol_displacement", "10")
#        self.rem.add("geom_opt_tol_energy", "5")

    def jobSopt(self):
        self.jobOpt()
        self.addSolvent()

    #Frequency Job
    def jobFreq(self):
        self.rem.jobtype("freq")
        self.rem.add("cpscf_nseg", "4")

    #standardized Q-Chem job naming scheme
    def genName(self):
        chmult = charge_tran[self.charge]+'m'+self.mult
        #self.name = self.path + self.name
        #print(self.name)
        name = self.path+'_'.join([self.name, self.jobtype, chmult])+'.in'
        return name

    @staticmethod
    def ReadChMult(chmult):
        ch, mult = chmult[:2], chmult[2:]
        ch = charge_tran[ch]
        mult = mult[1]
        return ch, mult

    def assembleJob(self):
        job = qc.inputfile()
        self.jobdict[self.jobtype]()
        job.add(self.rem)
        job.add(self.mol)
        for job_arr in self.job_arr_list:
            job.add(job_arr)
        return job

    def writeQcIn(self,qc_infile=None):
        job = self.assembleJob()
        print("QC infile: ", qc_infile)
        if qc_infile == None:
            qc_infile = self.genName()
        job.write(qc_infile)

def WriteJob(infile, args, indir='.'):
    if args.ParseName == True:
        aNmN = infile.split('.')[0].split('_')[-1]
        read_c, read_m = QcIn.ReadChMult(aNmN)
    else:
        read_c, read_m = args.c, args.m

    if indir == '.':
        readfile = infile
    #    qcw = QcIn(infile,read_c,read_m,args.j,basis=args.basis,method=args.method, \
    #               nametrunc=args.name)
    else:
        readfile = indir+'/'+infile
    #    qcw = QcIn(indir+'/'+infile,read_c,read_m,args.j,basis=args.basis,method=args.method, \
    #               nametrunc=args.name)
    qcw = QcIn(readfile,read_c,read_m,args.j,basis=args.basis,method=args.method, \
               nametrunc=args.name, infodump=args.infodump)
    qcw.writeQcIn()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='Q-Chem Output File or Parent Directory', type=str)
    parser.add_argument('-j', help='Job Type', type=str,default='sp')
    parser.add_argument('-c', help='Charge', type=str,default='0')
    parser.add_argument('-m', help='Multiplicity', type=str,default='1')
    parser.add_argument('-basis', help='Basis Set', type=str,default='6-31+g*')
    parser.add_argument('-method', help='EST Method', type=str,default='tpssh')
    parser.add_argument('-name', help='Truncate name at first _', type=int,default=1)
    parser.add_argument('-ParseName', help='Parse name for charge, spin mult', type=int,default=0)
    parser.add_argument('-infodump', help='add Molden format, save checkpoint', type=int, default=0)
    args = parser.parse_args()

    f = args.f
    if os.path.isfile(f):
        WriteJob(f, args)
#        if args.ParseName == True:
#            aNmN = f.split('.')[0].split('_')[-1]
#            read_c, read_m = QcIn.ReadChMult(aNmN)
#        qcw = QcIn(f,args.c,args.m,args.j,basis=args.basis,method=args.method,nametrunc=args.name)
#        qcw.writeQcIn()

    elif os.path.isdir(f):
        print("Directory found")
        for qcoutfile in os.listdir(f):
            if qcoutfile.split('.')[-1] in ('xyz', 'out'):
                print(qcoutfile)
                WriteJob(qcoutfile, args, f)
#                if args.ParseName == True:
#                    aNmN = qcoutfile.split('.')[0].split('_')[-1]
#                    read_c, read_m = QcIn.ReadChMult(aNmN)
#                    qcw = QcIn(f+'/'+qcoutfile,read_c,read_m,args.j,basis=args.basis,method=args.method, \
#                               nametrunc=args.name)
#                else:
#                    qcw = QcIn(f+'/'+qcoutfile,args.c,args.m,args.j,basis=args.basis,method=args.method, \
#                               nametrunc=args.name)
#                qcw.writeQcIn()

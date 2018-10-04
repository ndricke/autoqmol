#!/usr/bin/env python2.7
import ChemData as CD
import sys
import numpy as np

#Single point/general job parse keys
finish_key ='Thank you very much for using Q-Chem.  Have a nice day.'
sp_key ='Total energy in the final basis set'
coord_key ='Standard Nuclear Orientation (Angstroms)'

#Freq job keys
freq_key ='Frequency:   '
entropy_key ='Total Entropy:'
enthalpy_key ='Total Enthalpy:'
zpe_key = 'Zero point vibrational energy:'
vibH_key ='Vibrational Enthalpy:'

#PCM job keys
gelec_key ='G_electrostatic'
solG_key ='Total Free Energy (H0 + V/2 + non-elec)'
solEtot_key ='Total energy in the final basis set'
job_key = 'jobtype'
chmul_key = '$molecule'
Mulliken_key = 'Ground-State Mulliken Net Atomic Charges'
Chelpg_key = 'Ground-State ChElPG Net Atomic Charges'
optgeo_key = '**  OPTIMIZATION CONVERGED  **'
orbital_key = "Orbital Energies (a.u.)"

#Parsing class for Q-Chem output files
class Qdata(object):
  def __init__(self):
    self.chelpg_flag = False
    #There was a point where I didn't include all of these, but hey, why not? It runs fast enough
    self.parse_base = {finish_key:self.isComplete, \
                       chmul_key:self.chargeMult, \
                       sp_key:self.spEnergy, \
                       coord_key:self.coordGrab, \
                       entropy_key:self.entropy,\
                       enthalpy_key:self.enthalpy,\
                       zpe_key:self.zeroPointEnergy,\
                       vibH_key:self.enthalpyVib, \
                       gelec_key:self.solvationG,
                       solG_key:self.solvEnergy, \
                       job_key:self.jobtype, \
                       Mulliken_key : self.Mulliken, \
                       Chelpg_key : self.Chelpg, \
                       optgeo_key : self.optCoordGrab, \
                       #orbital_key : self.orbitals, \ #Currently breaks parsing for all things after the orbitals
                       }

#main workhorse; iterate through file and grab stuff when it hits a key phrase
  def qParse(self, qchem_outfile):
    with open(qchem_outfile, 'r') as f:
      for line in f:
        #if match with parse list, call function
        
        self.trash_key_bin = [] #Any keys to be removed will be removed after iterating over the dictionary
        for item in self.parse_base.keys():
          if item in line:
            func_point = self.parse_base[item]
            func_point(f,line)
            
        for trash_key in self.trash_key_bin:
            del self.parse_base[trash_key]
            print("Function parse routine complete. Removing key: "+trash_key)

  def readFile(self, filename):
    self.filename = filename.split('/')[-1]
    spl_filename = filename.split('.')
    self.ftype = spl_filename[-1]
    if self.ftype == 'xyz': #Just get the coordinates from the xyz file
      with open(filename,'r') as f:
        coord_list = f.read().splitlines(True)[2:]
        coord_list = [line.strip('\n') for line in coord_list]
        self.coordArr(coord_list)
    elif self.ftype == 'out': #Q-Chem output file
      self.qParse(filename)

  def default(self):
    self.basis = '6-31+g*'
    self.method = 'b3lyp'
    self.mult = 1
    self.charge = 0
    self.job = 'opt'
    #self.solvation = 'pcm'
    self.solvation = 'None'

  def listCoord(self):
    coord_list = []
    for i,atom in enumerate(self.atoms):
      coord_list.append(' '.join([atom] + [str(co) for co in list(self.coord[i,:])]))
    return coord_list

  def isComplete(self, infile = None,line=None):
    self.finish = True

  def chargeMult(self,infile,line):
    chmul = next(infile).strip('\n').split()
    self.charge = int(chmul[0])
    self.mult = int(chmul[1])
    #del self.parse_base[chmul_key]
    self.trash_key_bin.append(chmul_key)


  def Mulliken(self, infile, line):
    key_end = "--------------"
    atom_list, charge_list, spin_list = CD.ChargeGrab(key_end, infile)
    self.mulliken = charge_list
    self.spin = spin_list

  def Chelpg(self, infile, line):
    if self.chelpg_flag == False:
        key_end = "--------------"
        atom_list, charge_list, spin_list = CD.ChargeGrab(key_end, infile)
        self.chelpg = charge_list
        self.chelpg_flag = True

  def orbitals(self, infile, line):
    alpha_key = "Alpha MOs"
    beta_key = "Beta MOs"
#    for line in infile:
#      if alpha_key in line:
#        self.alpha_occ, self.alpha_vir = CD.moGrab(infile, line)
#      if beta_key in line:
#        self.beta_occ, self.beta_vir = CD.moGrab(infile, line)

  def spEnergy(self, infile, line):
    spline = line.split()
    self.E = float(spline[-1])

  def solvEnergy(self,infile,line):
    spline = line.split()
    self.E = float(spline[-2])

  def entropy(self, infile, line):
    spline = line.split()
    self.S = float(spline[2])/1000.0
    #del self.parse_base[entropy_key]

  def enthalpy(self, infile, line):
    spline = line.split()
    self.H = float(spline[2])
    #del self.parse_base[enthalpy_key]

  def enthalpyVib(self, infile, line):
    spline = line.split()
    self.Hvib = float(spline[2])
    #del self.parse_base[vibH_key]

  def frequencies(self, infile, line):
    spline = line.split()
    [self.freq.append(float(frq)) for frq in spline[1:]]

  def zeroPointEnergy(self, infile, line):
    spline = line.split()
    self.Hzpe = float(spline[-2])
    #del self.parse_base[zpe_key]

  def solvationG(self, infile, line):
    spline = line.split()
    self.Gpcm = float(spline[-2])
    self.solvation = 'pcm'

  def coordGrab(self, infile, line):
    key_end = "--------------"
    next(infile) ; next(infile) #Skipping junk. 2nd line contains ------
    grab_list = []
    for line in infile:
      if key_end in line:
        break
      grab_list.append(line.strip())
    self.atoms, self.coord = self.coordArr(grab_list)
  
  def jobtype(self, infile, line):
    spline = line.split()
    self.job = spline[-1]

  def optCoordGrab(self, infile, line):
    key_end = "Z-matrix Print:"
    next(infile); next(infile); next(infile); next(infile);
    grab_list = []
    for line in infile:
        if key_end in line:
            break
        grab_list.append(line.strip())
    print(grab_list)
    grab_list.pop(-1)

    self.opt_atoms, self.opt_coord = self.coordArr(grab_list)


  def coordArr(self,coord_list):
    if self.ftype == 'xyz': split_at = 1
    elif self.ftype == 'out': split_at = 2
    coord = np.zeros((len(coord_list),3))
    atoms = []
    for i, line in enumerate(coord_list):
      spline = line.split()
      atoms.append(spline[split_at-1])
      coord[i,:] = [float(j) for j in spline[split_at:]]
    return atoms, coord

  def dist(self,atom_ind1,atom_ind2):
    return np.linalg.norm(self.coord[atom_ind1,:] - self.coord[atom_ind2,:])

  @staticmethod
  def writeXYZ(filename, atom_list, coord_arr, comment=" "):
    with open(filename, 'w') as f:
      f.write(str(len(atom_list))+'\n')
      f.write(comment+'\n')
      for i, item in enumerate(atom_list):
        f.write(item + '   ' + '   '.join([str(num) for num in coord_arr[i,:]])+'\n')

    

    

#short testing module for reading a frequency calculation
if __name__ == "__main__":
  qdata = Qdata()
  qdata.readFile(sys.argv[1])

#  print(qdata.Hvib)
#  print(qdata.alpha_occ)
#  print(qdata.beta_occ)
#  print(len(qdata.alpha_occ))
#  print(len(qdata.beta_occ))

#  print(qdata.atoms)
#  print(qdata.chelpg)
  
#  Fe_charge = qdata.chelpg[qdata.atoms.index('Fe')]
#  print(Fe_charge)
#  print(qdata.spin)  

  print(qdata.opt_coord)
  qdata.writeXYZ('test.xyz', qdata.opt_atoms, qdata.opt_coord)





























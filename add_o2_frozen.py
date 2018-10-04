# input: directory containing bare catalyst .mol and/or .xyz files
# creates .xyz files for corresponding O2-bound structures
# catalysts with one metal atom only

import os, sys
import molSimplify
from molSimplify.Scripts import structgen
from molSimplify.Classes import mol3D, atom3D
import numpy as np
import math


class AddO2(object):

    def __init__(self, infile, active_sites=None):
        # also, molecule no longer explicitly optimized in this function
        self.infile = infile
        infile_type = self.infile.split('.')[-1]
        print(self.infile)

        self.molecule = mol3D.mol3D()
        self.molecule.OBMol = self.molecule.getOBMol(infile, infile_type, ffclean = False)
        self.molecule.convert2mol3D()

        if active_sites == None:
            self.active_sites = self.standardActiveSites()
        else:
            self.active_sites = active_sites

        self.bond_length_cutoffs = {'C': 1.953, 'H': 1.443, 'I': 2.543, 'Cl': 2.154,
                                   'B': 2.027, 'N': 1.834, 'F': 1.707, 'Br': 2.423,
                                   'P': 2.297, 'S': 2.198, 'O': 1.810} # X-O (angstroms)
        self.dtheta = 5 # change this if you want to

    def standardActiveSites(self):
        active_site_dict = {"mepyr": [14,16], "tetrid": [16], "tetry": [18,26]}
        metal = self.molecule.findMetal()
        if len(metal) >= 1:
            self.catO_bond_length = 1.8
            self.bond_angle = 120
            self.reposition_O2 = False
            return metal
        else:
            for key in active_site_dict.keys():
                if key in self.infile:
                    self.catO_bond_length = 1.55
                    self.bond_angle = 111
                    self.reposition_O2 = True
                    return active_site_dict[key]

    @staticmethod
    def scale_vector(vec, new_mag):
        squared = 0.0
        new_mag = float(new_mag)
        scale = new_mag / np.linalg.norm(vec)
        return [c * scale for c in vec]

    @staticmethod
    def vec3_perp(vec, theta):
        # janky orthogonal vector generator; 3D
        # theta = 0 (degrees) corresponds to i-hat direction. this is arbitrary
        if len(vec) != 3 or vec == [0, 0, 0]:
            print("Not a suitable vector.")
            return None
        v_theta = [math.cos(math.radians(theta)), math.sin(math.radians(theta)), 0]
        cross = np.cross(vec, v_theta)
        if [c for c in cross] == [0, 0, 0]: # could happen if vec and v_theta are parallel
            cross = np.cross(vec, v_theta + [0, 0, 1])
        return AddO2.scale_vector(cross, np.linalg.norm(vec)) # retains original magnitude

    def bind_O2(self, site_index, catO_bond_length=1.8, bond_angle=120, OO_bond_length=1.3, reposition_O2=False):
        # binds O2 perpendicular to catalyst plane, then makes .xyz file
            # does not modify input molecule
        # bond lengths in angstroms, angle in degrees. Variables set to parameters for Fe
        molecule_copy = mol3D.mol3D()
        molecule_copy.copymol3D(self.molecule)

        connection_list = molecule_copy.getBondedAtomsSmart(site_index, oct = False)
        #print("Active Site: %s       Neighbors: %s" %(str(site_index), str(connection_list))) # for debugging
        try:
            v1 = molecule_copy.atoms[connection_list[0]].distancev(molecule_copy.atoms[connection_list[1]])
            v2 = molecule_copy.atoms[connection_list[0]].distancev(molecule_copy.atoms[connection_list[2]])
            v = np.cross(v1, v2)
        except:
            print("Error finding plane of catalyst.")
            return None

        v_O1 = self.scale_vector(v, catO_bond_length)
        v_O2 = self.scale_vector(v, catO_bond_length + OO_bond_length * math.cos(math.pi - math.radians(bond_angle)))
        theta_0 = 0 # arbitrarily chosen
        v_perp = self.vec3_perp(self.scale_vector(v, OO_bond_length * math.sin(math.pi - math.radians(bond_angle))), theta_0)
        if v_perp == None: # shouldn't happen
            return None
        site_coords = molecule_copy.atoms[site_index].coords()

        molecule_copy.addAtom(mol3D.atom3D(Sym = 'O', xyz = [site_coords[c] + v_O1[c] for c in range(len(site_coords))]))
        molecule_copy.addAtom(mol3D.atom3D(Sym = 'O', xyz = [site_coords[c] + v_O2[c] + v_perp[c] for c in range(len(site_coords))]))
        #structgen.ffopt('MMFF94', molecule_copy, [], 1, [], False, [], 200, False) # keeps failing

        if reposition_O2: # workaround for ffopt() failing to set up; nonmetal catalysts
            theta = theta_0 + self.dtheta
            position_OK = False
            while not position_OK:
                position_OK = True
                for atom in molecule_copy.getAtoms():
                    if atom == molecule_copy.getAtoms()[-1] or atom == molecule_copy.getAtoms()[-2]:
                        continue
                    if atom.symbol() in self.bond_length_cutoffs:
                        cutoff = self.bond_length_cutoffs[atom.symbol()]
                    else: # won't happen if bond_length_cutoffs is kept up to date
                        cutoff = self.bond_length_cutoffs[max(self.bond_length_cutoffs, key = lambda key: self.bond_length_cutoffs[key])]
                        print("Using %f A as the bond length cutoff for %s" %(cutoff, atom.symbol()))
                    if molecule_copy.getAtoms()[-1].distance(atom) < cutoff:
                        position_OK = False
                        break
                if not position_OK:
                    if theta >= 360:
                        print("%s may have geometry issues after O2-binding." %infile)
                        break
                    v_perp = self.vec3_perp(self.scale_vector(v, OO_bond_length * math.sin(math.pi - math.radians(bond_angle))), theta)
                    if v_perp == None: # shouldn't happen
                        return None
                    molecule_copy.getAtoms()[-1].setcoords([site_coords[c] + v_O2[c] + v_perp[c] for c in range(len(site_coords))])
                    theta += self.dtheta

        file_name = self.infile.split('.')[0] + "O2"
        if len(self.active_sites) > 1: #only add active site in name if more than 1
            file_name += "-" + str(site_index) # still zero-indexed
        molecule_copy.writexyz(file_name)

    def run(self):
        for site in self.active_sites:
            self.bind_O2(site, catO_bond_length=self.catO_bond_length, bond_angle=self.bond_angle, reposition_O2=self.reposition_O2)

if __name__ == "__main__":
    infile_dir = sys.argv[1]
    if os.path.isdir(infile_dir):
        file_list = os.listdir(infile_dir)
        for infile in file_list:
            molecule_O2 = AddO2(infile_dir + '/' + infile)
            molecule_O2.run()
    else:
        molecule_O2 = AddO2(infile)
        molecule_O2.run()

#!/usr/bin/env python2.7

import pyQChem as qc
import sys
import os
from qcIn import QcIn


"""
This is supposed to make it easy to build consecutive jobs
To do that, we want some flexibility about the types of jobs that can be run consecutively
*assuming charge, mult,
"""

class QcMultIn(object):

    def __init__(self, infile, charge, mult, joblist, method='tpssh', nametrunc=True, basis='6-31+g*'):
        self.charge_tran = {'a2':'-2','a1':'-1','a0':'0','c1':'1','c2':'2','c3':'3','c4':'4'}
        revd=dict([reversed(i) for i in self.charge_tran.items()])
        self.charge_tran.update(revd)

        splinfile = infile.split('/')
        if len(splinfile) == 1:
            self.path = ''
        else:
            self.path = '/'.join(splinfile[:-1])+'/'
        self.name = (splinfile[-1].split('.'))[0]
        if nametrunc == True:
            self.name = self.name.split('_')[0]
        self.infile = infile
        self.charge = charge
        self.mult = mult
        print(self.charge, self.mult)

        self.joblist = joblist #joblist is a list of Q-Chem jobtypes
        self.method = method
        self.jobconcat = ''.join(joblist)
        self.job = qc.multifile()
        self.nametrunc = nametrunc
        self.basis = basis

    def assembleMulti(self):
        qcw = QcIn(self.infile,self.charge,self.mult, self.joblist.pop(0), method=self.method, nametrunc=self.nametrunc, basis=self.basis)
        self.job.add(qcw.assembleJob())
        for jobitem in self.joblist:
            qcw = QcIn(self.infile,self.charge,self.mult, jobitem, method=self.method, geom_read=True, nametrunc=self.nametrunc, basis=self.basis)
            self.job.add(qcw.assembleJob())

    def writeMulti(self):
        self.assembleMulti()
        #qc_infile = self.genName()
        qc_infile = self.path + self.genName()
        self.job.write(qc_infile)

    #standardized Q-Chem job naming scheme for multi-jobs
    def genName(self):
        chmult = self.charge_tran[self.charge]+'m'+self.mult
        name = '_'.join([self.name, self.jobconcat, chmult])+'.in'
        return name

def chargeTran():
        charge_tran = {'a2':'-2','a1':'-1','a0':'0','c1':'1','c2':'2','c3':'3','c4':'4'}
        revd=dict([reversed(i) for i in charge_tran.items()])
        charge_tran.update(revd)
        return charge_tran

def chmultName(qcoutname):
        charge_tran = chargeTran()
        name_split = qcoutname.split('.')[0].split('_')
        chmult = name_split[-1].split('m')
        assert len(chmult) == 2
        c = charge_tran[chmult[0]]
        m = chmult[-1]
        return c, m

if __name__ == "__main__":
    import argparse
    import copy

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='Q-Chem Output File or Parent Directory', type=str)
    parser.add_argument('-j', help='Job Type', type=str,default='optsp')
    parser.add_argument('-c', help='Charge', type=str,default='0')
    parser.add_argument('-m', help='Multiplicity', type=str,default='1')
    parser.add_argument('-basis', help='Basis Set', type=str,default='6-31+g*')
    parser.add_argument('-method', help='EST Method', type=str,default='b3lyp')
    parser.add_argument('-nameparse', help='EST Method', type=bool,default=False)
    parser.add_argument('-nametrunc', help='Truncate name at first _', type=int,default=1)
    args = parser.parse_args()

    multijob_dict = {'optsp':['opt','sp'], 'cdftsp':['cdft','sp'], 'optcdftsp':['opt','cdft','sp']}
    job = multijob_dict[args.j]

    c, m = args.c, args.m
    f = args.f
    if os.path.isfile(f):
        if args.nameparse:
            c,m = chmultName(f)
        qcw = QcMultIn(f,c,m, job, method=args.method, nametrunc=args.nametrunc, basis=args.basis)
        qcw.writeMulti()

    elif os.path.isdir(f):
        for qcoutfile in os.listdir(f):
#            if qcoutfile.split('.')[-1] == 'out': ##PyQchem has some trouble with certain Q-chem output files, so xyz is safer
            if qcoutfile.split('.')[-1] == 'xyz':
                if args.nameparse:
                    c,m = chmultName(qcoutfile)
                qcw = QcMultIn(f+'/'+qcoutfile, c, m, copy.copy(job), method=args.method, nametrunc=args.nametrunc, basis=args.basis)
                qcw.writeMulti()

#!/usr/bin/env python
#coding=utf-8
# 添加了iterator()方法用于遍历当前目录下的文件夹中的XDATCAR文件(6-4-2016)

from matplotlib.font_manager import FontProperties
import gc
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import os
import math
import copy
from progressbar import *

class Extract():
    # ensemble: ensemble type of md. e.g. nvt, npt
    def __init__(self, ensemble="npt"):
        self.ensemble = ensemble

    # identifier: indicate the start of a structure to distinguish these different structures
    def read(self, identifier):
        os.system("grep 'energy  without entropy=' OUTCAR > e.txt")
        
        infile = open("XDATCAR")
        infile2 = open("e.txt")

        #outfile = open("input.data", "w")
        structures = []
        energys = []
        
        string = infile.readline()
        while(string):
            if (self.ensemble == "npt"):
                if (string.startswith(identifier)):
                    scale = float(infile.readline()) # scale of lattice parameter
                    # lattice parameter
                    a = []
                    for i in xrange(0, 3):
                        string = infile.readline()
                        tmp = [float(s0) for s0 in string.split()]
                        a.append(tmp)

                    a = np.array(a)*scale

                    # type of element
                    string = infile.readline()
                    element = np.array(string.split())
                    string = infile.readline()
                    elementNumber = np.array([int (s0) for s0 in string.split()])

                    # read atom coordinate
                    if (string.startswith("Direct configuration=")):
                        atoms = self.atomCoordinate(infile, a, element, elementNumber)
                        tmp = np.vstack([a, atoms])
                        tmp = list(tmp) # array of numpy to list of python
                        structures.append(tmp)
                        energys.append(float(infile2.readline().splie()[3]))
                        
            elif (self.ensemble == "nvt"):
                if (string.startswith(identifier)):
                    scale = float(infile.readline()) # scale of lattice parameter
                    # lattice parameter
                    a = []
                    for i in xrange(0, 3):
                        string = infile.readline()
                        tmp = [float(s0) for s0 in string.split()]
                        a.append(tmp)

                    a = np.array(a)*scale
                    
                    # type of element
                    string = infile.readline()
                    element = np.array(string.split())
                    string = infile.readline()
                    elementNumber = np.array([int (s0) for s0 in string.split()])
                    
                # read atom coordinate
                if (string.startswith("Direct configuration=")):
                    atoms = self.atomCoordinate(infile, a, element, elementNumber)
                    tmp = np.vstack((a, atoms))
                    tmp = list(tmp) # array of numpy to list of python
                    structures.append(tmp)
                    energys.append(float(infile2.readline().split()[3]))

            string = infile.readline()
        structures = np.array(structures)
        energys = np.array(energys)

        return structures, energys

    # read atom coordinate of a structure
    def atomCoordinate(self, infile, a, element, elementNumber):

        atoms = [] # coordinate of atoms
        for i in xrange(0, np.sum(elementNumber)):
            string = infile.readline()
            tmp = np.array([float(s0) for s0 in string.split()])
            tmp = self.direct2cartesian2(a, tmp) # direct to cartesian
            atoms.append(tmp)
        atoms = np.array(atoms)
        
        return atoms

    # 坐标转换
    # 分数坐标转换为直角坐标
    def direct2cartesian(self, a, coordinate):
        a = copy.deepcopy(a)
        coordinate = copy.deepcopy(coordinate)

        for i in xrange(0, coordinate.shape[0]): # atom
            tmp = 0
            for j in xrange(0, coordinate.shape[1]): # dirction
                tmp += coordinate[i][j]*a[j]
            coordinate[i] = tmp
        return coordinate

    # only for one coordiante
    def direct2cartesian2(self, a, coordinate):
        a = copy.deepcopy(a)
        coordinate = copy.deepcopy(coordinate)

        tmp = 0
        for i in xrange(0, coordinate.shape[0]): # dirction
            tmp += coordinate[i]*a[i]
        coordinate = tmp
        return coordinate

    # iterator is used to extact all directoris
    # files: [start, end]
    # filesLength: length of file name
    # identifier: indicate the start of a structure to distinguish these different structures
    def iterate(self, files, fileLength, identifier):
        path = os.getcwd()

        # 进度条
        pbar = ProgressBar(widgets=['Read: ', Percentage(), ' ', Bar(),
                                    ' ', ETA()], maxval=files[1]).start()
        
        structures = []
        energys = []
        for ii in xrange(files[0], files[1]):
            dir = (fileLength - len(str(ii))) * "0" + str(ii)
       
            if os.path.isdir(path+"/"+dir):
                os.chdir(path+'/'+dir)
                try:
                    tempS, tempE = self.read(identifier)
                except IOError:
                    print "XDATCAR is inexistent!", ii
                    os.chdir(path)
                    break

                if (structures == []):
                    structures = tempS
                    energys = tempE
                else:
                    structures = np.concatenate((structures, tempS), axis=0)
                    energys = np.concatenate((energys, tempE), axis=0)
                #print os.getcwd()
                os.chdir(path)
            # 进度条
            #time.sleep(0.01)
            pbar.update(ii)
            

            
        pbar.finish()
        return structures, energys
        
    # write structure and its energy to file
    def outStructure(self, outfile, structure, energy):
        outfile.write("Next Structure\n")

        for i in xrange(0, structure.shape[0]):
            if (i < 3):
                outfile.write("{0:.6f}  {1:.6f}  {2:.6f}\n".format(structure[i][0], structure[i][1], structure[i][2]))
            else:
                outfile.write("{0:.12f}  {1:.12f}  {2:.12f}\n".format(structure[i][0], structure[i][1], structure[i][2]))
        # energy
        outfile.write("{0:.8f}\n".format(energy))
        
    # output structures
    # srange: the range of structure, [start, end ,step]
    def output(self, structures, energys, srange=None):
        outfile = open("input.data", "w")
        if (srange == None):
            for i in xrange(0, structures.shape[0]):
                self.outStructure(outfile, structures[i], energys[i])
        else:
            for i in xrange(srange[0], srange[1], srange[2]):
                self.outStructure(outfile, structures[i], energys[i])
                
        outfile.close()
        
# --------------------

e = Extract(ensemble="nvt")
# md
#structures, energys = e.read("Te")
#print structures.shape, energys.shape
#e. output(structures, energys, [0, energys.shape[0], 2])

# thirdorder_vasp
structures, energys = e.iterate([1,1013], 4, "SiO2")
print structures.shape, energys.shape
print os.getcwd()
e. output(structures, energys, [0, energys.shape[0], 1])

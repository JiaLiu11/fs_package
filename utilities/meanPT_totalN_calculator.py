#!/usr/bin/env python 

#   file: meanPT_totalN_calculator.py
#
#   Author:   Jia Liu    <liu.2053@osu.edu>
#   History:
#   Jul 16, 2014     First version.


import sys, shutil
from os import path, stat, getcwd, stat, rename, remove, makedirs
import numpy as np
from subprocess import call

# run parameters
event_list = range(1,41)
matchingTime_list = np.linspace(0.4,1.8,8)

# folders locations
rootDir = path.abspath('..')  #root directory for current node
database_location = path.join(path.abspath('../../..'),'dataBase') # 
fs_location = path.join(rootDir, 'fs')
fs_particle_location = path.join(rootDir, 'fs','fs_particle_data')
is_location = path.join(rootDir, 'iS')
table_location = path.join(rootDir, 'tables')
node_name = rootDir.split('/')[-1]
sfactor_list = np.loadtxt(path.join(table_location,'sfactor_log.dat'))


def readParticleNum(fileName, targetFolder):
    """
    read particle number from file.
    return: total particle number
    """
    fileName = path.join(targetFolder, fileName)
    try:
        particle_data = np.loadtxt(fileName)
    except:
        print 'dEcounters: readParticleNum error!'
        print 'Cannot read in particle file: ' + fileName
        sys.exit(-1)
    return particle_data[0, 1]


def calculateParticleMeanPT(fileName, pT_tbl, pTweight_tbl, targetFolder):
    """
    calculate the mean pT of a thermal particle or all charged particle.
    return: (mean pT, total number) of a thermal particle, 
    """
    # find out if hydro ran
    decdat2_file = path.join(targetFolder, 'decdat2.dat')
    if(stat(decdat2_file).st_size==0):  # empty decdat2 file, hydro did not run
        return (0, 0)   
    # read in particle dN/(dyptdpt2\pi) table
    fileName = path.join(targetFolder, fileName)
    try:
        particle_data = np.loadtxt(fileName)
    except:
        print 'dEcounters: calculateParticleMeanPT() error!'
        print 'Cannot read in particle file: ' + fileName
        sys.exit(-1)
    dndyptdpt_2pi_tbl = particle_data[:, 2]
    # get particle pt
    total_pt = (dndyptdpt_2pi_tbl * pT_tbl * pT_tbl * pTweight_tbl * 2.0 * np.pi).sum()
    return total_pt



def meanPTCalculatorShell():
    # prepare integration table
    pt_tbl_file = path.join(rootDir, 'iS/tables', 'pT_gauss_table.dat')
    pt_tbl = np.loadtxt(pt_tbl_file)
    meanPT_file = path.join(database_location, node_name, 'meanPT_totalN.dat')
    meanPT_log = open(meanPT_file, 'w')

    for event_num in event_list:
        for tau_s in matchingTime_list:
            # get photon pT
            is_data_folder = path.join(database_location, node_name, 'event_%d'%event_num,'%g'%tau_s)
            allch_integrated_filename = 'Charged_ptcut015_10_eta_integrated_vndata.dat'
            allch_filename = 'Charged_ptcut015_10_eta_vndata.dat'
            allch_num = readParticleNum(allch_integrated_filename, is_data_folder)
            allch_pt = calculateParticleMeanPT(allch_filename,pt_tbl[:,0], pt_tbl[:,1],is_data_folder)
            # consider the degenercy of gluon and photon
            print "%8.2f \t %10.6e \t %10.6e \t %10.6e "%(tau_s, allch_pt, allch_num, allch_pt/allch_num)
            meanPT_log.write("%8.2f \t %10.6e \t %10.6e \t %10.6e "%(tau_s, allch_pt, allch_num, allch_pt/allch_num))
            meanPT_log.flush()
    meanPT_log.close()

if __name__ == "__main__":
    meanPTCalculatorShell()

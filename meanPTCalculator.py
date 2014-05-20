#!/usr/bin/env python 

#   file: meanPTCalculator.py
#
#   Author:   Jia Liu    <liu.2053@osu.edu>
#   History:
#   May 19, 2014     first version

import sys, shutil
from os import path, stat, getcwd
import numpy as np
from subprocess import call

rootDir = getcwd()
fs_location = path.join(rootDir, 'fs')
fs_particle_location = path.join(rootDir, 'fs_particle')
table_location = path.join(rootDir, 'tables')
sfactor_list = np.loadtxt(path.join(table_location,'sfactor_log.dat'))

def generatePartonNumFromLm(datafile, event_num, tau0, taumin, taumax, dtau):
    """
    Use free-streaming and matching code to process parton pT^1 file and dump the dNdydphi file.
    Return: success code or fail code.
    """
    # get data file and copy it to fs_folder
    fs_init_folder = path.join(fs_particle_location, 'data', 'events')
    fs_init_file = path.join(fs_init_folder, 'sd_event_%d_block.dat' % event_num)
    shutil.copy2(datafile, fs_init_file)
    # run free-streaming and matching
    lm_cmd = fs_particle_location + "/./lm.e " + 'event_mode=' + str(event_num) \
        + ' tau0=' + str(tau0) + ' taumin=' + str(taumin) + ' taumax=' + str(taumax) \
        + ' dtau=' + str(dtau)
    lm_retcode = call(lm_cmd, shell=True, cwd=fs_particle_location)
    if lm_retcode == 0:
        pass
    else:
        print 'generatePartonNum: Freestreaming and Landau Matching stops unexpectly!'
        sys.exit(-1)
    # get particle number
    return


def calculatePartonMeanPT(event_num, tau_s, edinit_file, sfactor, Edec, dxdy=0.01):
    """
    calculate mean pT of un-thermalized partons.
    Return: (mean pT, total number) of massless particles.
    """
    # identify necessary files
    dptd2rdphi_file = path.join(fs_location, 'data', 'result', 'event_%d'%event_num, '%g'%tau_s, 'dEd2rdphip_kln.dat')
    dnd2rdphi_file = path.join(fs_particle_location, 'data', 'result', 'event_%d'%event_num, '%g'%tau_s,\
        'dEd2rdphip_kln.dat')
    phipTbl = np.loadtxt(path.join(table_location,'phip_gauss_table.dat'))  #100 points Gaussian points table for phip integration.
    try:
        dptd2rdphi_data = np.loadtxt(dptd2rdphi_file)
        dnd2rdphi_data = np.loadtxt(dnd2rdphi_file)
        ed_data = np.loadtxt(edinit_file)
    except:
        print 'meanPTCalculator: calculatePartonMeanPT() error!'
        print 'No parton data file!'
        sys.exit(-1)
    phipGaussWeight = phipTbl[:, 1]
    ed_data = sfactor * ed_data
    ed_criteria = ed_data < Edec
    ed_criteria_rowNum = np.reshape(ed_criteria, (ed_data.size))  #reshape is done row-wise
    # non-zero rows ed<Edec
    # get dptdy
    dptd2rdphi_data = dptd2rdphi_data * sfactor
    dptd2rdphi_data_outside = dptd2rdphi_data[ed_criteria_rowNum, :]
    dptdydphip = dptd2rdphi_data_outside.sum(axis=0) * dxdy  # summing along each columnn
    dptdy = (dptdydphip * phipGaussWeight).sum()
    
    # get dndy
    dnd2rdphi_data = dnd2rdphi_data * sfactor
    dnd2rdphi_data_outside = dnd2rdphi_data[ed_criteria_rowNum, :]
    dndydphip = dnd2rdphi_data_outside.sum(axis=0) * dxdy  # summing along each columnn
    dndy = (dndydphip * phipGaussWeight).sum()
    
    return (dptdy, dndy)


def calculateParticleMeanPT(particle_idx, dpT_tbl, pTweight_tbl, targetFolder):
    """
    calculate the mean pT of a thermal particle.
    return: (mean pT, total number) of a thermal particle, 
    """
    # read in particle dN/(dyptdpt2\pi) table
    fileName = path.join(targetFolder, 'thermal_%d_vndata.dat' % particle_idx)
    try:
        particle_data = np.loadtxt(fileName)
    except:
        print 'dEcounters: calculateParticleMeanPT() error!'
        print 'Cannot read in particle file: ' + fileName
        sys.exit(-1)
    dndyptdpt_2pi_tbl = particle_data[:, 2]
    # get total particle number
    total_num = readParticleNum(particle_idx, targetFolder)
    # get particle pt
    total_pt = (dndyptdpt_2pi_tbl * dpT_tbl * dpT_tbl * pTweight_tbl * 2.0 * np.pi).sum()
    return (total_pt, total_num)


def readParticleNum(particle_idx, targetFolder):
    """
    read particle number from file.
    return: total particle number
    """
    fileName = path.join(targetFolder, 'thermal_%d_integrated_vndata.dat' % particle_idx)
    try:
        particle_data = np.loadtxt(fileName)
    except:
        print 'dEcounters: readParticleNum error!'
        print 'Cannot read in particle file: ' + fileName
        sys.exit(-1)
    return particle_data[0, 1]


def meanPTCalculatorShell():
    event_num = 99
    particle_file = path.join('superMC','data', 'sd_event_%d_block_particle.dat'%event_num)
    generatePartonNumFromLm(particle_file, event_num, 0.01, 1, 10, 1)
    print '     tau_s     total pT      total num      mean pT'
    for tau_s in range(1, 3):
        edinit_file = path.join('fs','data','result','event_%d'%event_num, '%g'%tau_s,'ed_profile_kln.dat')
        sfactor = (sfactor_list[sfactor_list[:,0]==tau_s])[0,1]
        totalpt, totalnum = calculatePartonMeanPT(event_num, tau_s, edinit_file, sfactor, 0.18)
        print "%8.2f \t %10.6e \t %10.6e \t %10.6e"%(tau_s, totalpt, totalnum, totalpt/totalnum)

if __name__ == "__main__":
    meanPTCalculatorShell()

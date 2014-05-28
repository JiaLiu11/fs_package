#!/usr/bin/env python 

#   file: meanPTCalculator.py
#
#   Author:   Jia Liu    <liu.2053@osu.edu>
#   History:
#   May 28, 2014     Get photon instead of pion pT from hydro, and make the weighting of 
#                   gluon:photon=1:1
#   May 19, 2014     first version

import sys, shutil
from os import path, stat, getcwd, stat, rename, remove
import numpy as np
from subprocess import call

rootDir = path.abspath('..')  #since this script is in fs_package/utilities/
fs_location = path.join(rootDir, 'fs')
fs_particle_location = path.join(rootDir, 'fs_particle')
is_location = path.join(rootDir, 'iS')
table_location = path.join(rootDir, 'tables')
matchingTime_list = np.linspace(1, 10, 10)
# some constants
photon_degeneracy = 2.0
gluon_degeneracy = 8.0

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
        + ' dtau=' + str(dtau) + ' >fs_runlog.dat'
    lm_retcode = call(lm_cmd, shell=True, cwd=fs_particle_location)
    if lm_retcode == 0:
        pass
    else:
        print 'generatePartonNum: Freestreaming and Landau Matching stops unexpectly!'
        sys.exit(-1)
    # get particle number
    return

def generatePartonEnergyFromLm(datafile, event_num, tau0, taumin, taumax, dtau):
    """
    Use free-streaming and matching code to process parton pT^2 file and dump the dNdydphi file.
    Return: success code or fail code.
    """
    # get data file and copy it to fs_folder
    fs_init_folder = path.join(fs_location, 'data', 'events')
    fs_init_file = path.join(fs_init_folder, 'sd_event_%d_block.dat' % event_num)
    shutil.copy2(datafile, fs_init_file)
    # run free-streaming and matching
    lm_cmd = fs_location + "/./lm.e " + 'event_mode=' + str(event_num) \
        + ' tau0=' + str(tau0) + ' taumin=' + str(taumin) + ' taumax=' + str(taumax) \
        + ' dtau=' + str(dtau) + ' >fs_runlog.dat'
    lm_retcode = call(lm_cmd, shell=True, cwd=fs_location)
    if lm_retcode == 0:
        pass
    else:
        print 'generatePartonNum: Freestreaming and Landau Matching stops unexpectly!'
        sys.exit(-1)
    # get particle number
    return

def calculatePartonMeanPT(event_num, tau_s, sfactor, Edec, dxdy=0.01):
    """
    calculate mean pT of un-thermalized partons.
    Return: (mean pT, total number) of massless particles.
    """
    # identify necessary files
    dptd2rdphi_file = path.join(fs_location, 'data', 'result', 'event_%d'%event_num, '%g'%tau_s, 'dEd2rdphip_kln.dat')
    dnd2rdphi_file = path.join(fs_particle_location, 'data', 'result', 'event_%d'%event_num, '%g'%tau_s,\
        'dEd2rdphip_kln.dat')
    edinit_file = path.join(fs_location, 'data', 'result', 'event_%d'%event_num, '%g'%tau_s, 'ed_profile_kln.dat')
    phipTbl = np.loadtxt(path.join(table_location,'phip_gauss_table.dat'))  #100 points Gaussian points table for phip integration.
    if (path.isfile(dptd2rdphi_file) and path.isfile(dnd2rdphi_file) and path.isfile(edinit_file)):
        dptd2rdphi_data = np.loadtxt(dptd2rdphi_file)
        dnd2rdphi_data = np.loadtxt(dnd2rdphi_file)
        ed_data = np.loadtxt(edinit_file)
    else:
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


################################################################################################
def runiSGetPhoton(data_location):
    """
    re-run iS only to find the photon pT spectra
    return: success code or fail code
    """
    # find out if hydro ran
    decdat2_file = path.join(data_location, 'decdat2.dat')
    if(stat(decdat2_file).st_size==0):  # empty decdat2 file, hydro did not run
        return 

    iS_EOS_location = path.join(is_location, 'EOS')
    # backup the choosen_particles table
    choosen_particles_table_original = path.join(iS_EOS_location, 'chosen_particles.dat')
    choosen_particles_table_backup = path.join(iS_EOS_location, 'chosen_particles.dat_backupi')
    if path.isfile(choosen_particles_table_original):
        rename(choosen_particles_table_original, choosen_particles_table_backup)
    # write only photon in the choosen particles file
    choosen_particles_table_photon = open(choosen_particles_table_original, 'w')
    choosen_particles_table_photon.write('  22\n')
    choosen_particles_table_photon.close()
    # move data files to iS
    decdat2_source = path.join(data_location, 'decdat2.dat')
    decdat2_target = path.join(is_location, 'results','decdat2.dat')
    shutil.copy2(decdat2_source, decdat2_target)
    surface_source = path.join(data_location, 'surface.dat')
    surface_target = path.join(is_location, 'results','surface.dat')
    shutil.copy2(surface_source, surface_target)
    decdatmu_source = path.join(data_location, 'decdat_mu.dat')
    decdatmu_target = path.join(is_location, 'results','decdat_mu.dat')
    shutil.copy2(decdatmu_source, decdatmu_target)
    # run iS
    iS_cmd = is_location + "/./iS.e >runlog.dat"
    retcode = call(iS_cmd, shell=True, cwd=is_location)
    if retcode ==0:
      pass
    else :
      print 'iS stops unexpectly!\n'
      sys.exit(-1)
    # resume the original choosen_particles file
    remove(choosen_particles_table_original)
    if path.isfile(choosen_particles_table_backup):
        rename(choosen_particles_table_backup, choosen_particles_table_original)
    # backup the photon file to database
    backup_cmd = 'cp results/thermal_22_*.dat '+ data_location
    call(backup_cmd, shell=True, cwd=is_location)
    return

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


def calculateParticleMeanPT(particle_idx, pT_tbl, pTweight_tbl, targetFolder):
    """
    calculate the mean pT of a thermal particle.
    return: (mean pT, total number) of a thermal particle, 
    """
    # find out if hydro ran
    decdat2_file = path.join(targetFolder, 'decdat2.dat')
    if(stat(decdat2_file).st_size==0):  # empty decdat2 file, hydro did not run
        return (0, 0)   
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
    total_pt = (dndyptdpt_2pi_tbl * pT_tbl * pT_tbl * pTweight_tbl * 2.0 * np.pi).sum()
    return (total_pt, total_num)



def getPhotonPT(dpT_tbl, pTweight_tbl, is_result_folder):
    """
    This is the shell of function calculateParticleMeanPT(). It gets mean pT of photon. 
    Input: tau_s, folder to the event 
    Return: total pT, total pion number
    """
    photon_MC_idx = 22 # Monte-Carlo number of photon
    total_pt = 0
    total_num = 0
    photon_folder = is_result_folder
    total_pt, total_num = calculateParticleMeanPT(photon_MC_idx, dpT_tbl, pTweight_tbl, photon_folder)
    return (total_pt, total_num)


def cleanfsResultFolders():
    """
    clean the fs results folder for this mean pt run.
    """
    fs_result_location = path.join(fs_location, 'data','result')
    fs_particle_result_location = path.join(fs_particle_location, 'data','result')
    shutil.rmtree(fs_result_location)
    shutil.rmtree(fs_particle_result_location)
    return


def meanPTCalculatorShell():
    # prepare integration table
    pt_tbl_file = path.join(rootDir, 'iS/tables', 'pT_gauss_table.dat')
    pt_tbl = np.loadtxt(pt_tbl_file)
    sfactor_list = np.loadtxt(path.join(table_location,'sfactor_log.dat'))

    event_num = 99
    MCevents_folder = path.join(rootDir,'superMC','data')
    particle_file = path.join(MCevents_folder, 'sd_event_%d_block_particle.dat'%event_num)
    energy_file   = path.join(MCevents_folder, 'sd_event_%d_block.dat'%event_num)
    # generate free-streamed profile again
    print 'Start to run free-streaming for pT order 1:'
    generatePartonNumFromLm(particle_file, event_num, 0.01, \
        matchingTime_list[0], matchingTime_list[-1], matchingTime_list[1]-matchingTime_list[0])
    print 'Start to run free-streaming for pT order 2:'
    generatePartonEnergyFromLm(energy_file, event_num, 0.01, \
        matchingTime_list[0], matchingTime_list[-1], matchingTime_list[1]-matchingTime_list[0])
    print 'fs code finished!'

    print '%    tau_s     total parton pT  total parton num total pion pT  total pion num    mean pT'
    for tau_s in matchingTime_list:
        # get parton pT
        sfactor = (sfactor_list[sfactor_list[:,0]==tau_s])[0,1]
        parton_totalpt, parton_totalnum = calculatePartonMeanPT(event_num, tau_s, sfactor, 0.18)
        # get photon pT
        is_data_folder = path.join(rootDir, 'localdataBase', 'event_%d'%event_num,'%g'%tau_s)
        runiSGetPhoton(is_data_folder)
        photon_pt, photon_totalnum = getPhotonPT(pt_tbl[:,0], pt_tbl[:,1],is_data_folder)
        # consider the degenercy of gluon and photon
        meanPT_all = (parton_totalpt/gluon_degeneracy+photon_pt/photon_degeneracy)\
        /(parton_totalnum/gluon_degeneracy+photon_totalnum/photon_degeneracy)
        print "%8.2f \t %10.6e \t %10.6e \t %10.6e \t %10.6e \t %10.6e"%(tau_s, parton_totalpt, parton_totalnum, \
            photon_pt, photon_totalnum, meanPT_all)
    # clean the temp files before leaving
    cleanfsResultFolders()

if __name__ == "__main__":
    meanPTCalculatorShell()

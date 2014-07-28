#!/usr/bin/env python 

#   file: runHydroOnly.py
#
#   Author:   Jia Liu    <liu.2053@osu.edu>
#
#   This python script only run hydro to update the avg_points.dat file.
#
#   Revision history:
#   Apr.24,2014  first version

import sys
from os import path, getcwd, listdir, rename, makedirs, remove, stat, popen
import shutil #contains command to move files
from subprocess import call
import numpy as np
import dEcounters

event_num_list = range(1, 41)
events_total = len(event_num_list)
matching_time_list = np.linspace(1, 10, 10)
sys_start_time = 0.01

angle_phi2 = 0.0
norm_factor = 1.0

rootDir = path.abspath('..')  #define the root folder

#define the data folder of superMC
superMCDirectory = path.join(rootDir, 'superMC')
superMCResultDirectory = path.join(superMCDirectory, 'data')

#define the data folder of fs
lmDirectory = path.join(rootDir, 'fs')  
lmInitDirectory = path.join(lmDirectory, 'data', 'events') 

#define the inital data folder of VISHNEW code
hydroDirectory = path.join(rootDir, 'VISHNew')
hydroInitialDirectory = path.join(hydroDirectory, 'Initial')
hydroResultDirectory = path.join(hydroDirectory, 'results')

# if run for multiple nodes
backupDirName = rootDir.split('/')[-1] # get the name of the current node
backupDir = path.join('..','..','..', 'dataBase', backupDirName)
tablesLocation = path.join(rootDir, 'tables')

sfactor_list = np.loadtxt(path.join(tablesLocation,'sfactor_log.dat')) #load in scale factor




def moveMCdata2Lm(source_dir, dest_dir):
    """
    move the initial profiles from superMC to fs data folder
    """
    cleanUpFolder(dest_dir)
    for init_files in listdir(source_dir):
        src_file = path.join(source_dir, init_files)
        dst_file = path.join(dest_dir, init_files)
        shutil.move(src_file, dst_file)
    print 'superMC profiles have been put into fs folder!' 

def runlm(event_num, tau0, taumin, taumax, dtau):
    """
    Run Landau matching code for a single event
    """
    lm_cmd = lmDirectory + "/./lm.e " + 'event_mode='+str(event_num) \
        + ' tau0='+str(tau0) +' taumin='+str(taumin) + ' taumax=' +str(taumax) \
        + ' dtau='+str(dtau)
    lm_retcode = call(lm_cmd, shell=True, cwd=lmDirectory)
    if lm_retcode ==0:
      print 'Freestreaming and Landau Matching completes!'
    else :
      print 'Freestreaming and Landau Matching stops unexpectly!'
      sys.exit(-1)




def moveLmData2Hydro(source_dir, dest_dir, tag):
    """
    move the initial conditions from landau matching to hydro
    """
    for lm_files in listdir(source_dir):
        new_file_name = lm_files.replace('_tauf_' + tag ,"")
        rename(path.join(source_dir,lm_files), 
               path.join(source_dir,new_file_name))
        #print lm_files, '->', new_file_name, 'rename ok'

    #move the renamed files to hydro
    for lm_files in listdir(source_dir):
        src_file = path.join(source_dir, lm_files)
        dst_file = path.join(dest_dir, lm_files)
        shutil.move(src_file, dst_file)
    #  print src_file, "->", dst_file
    print 'Initial conditions have been put into hydro folder\n'

def getSfactor(mth_time):
    """ 
    Get the rescale factor from the data file, which is read in at the 
    beginning of this script
    """
    # #find if the matching time is in the table
    # temp = '%.8f'%mth_time
    # mth_time = float(temp)
    # if (sfactor_list[:,0]==mth_time).any()==True:
    #     return (sfactor_list[sfactor_list[:,0]==mth_time])[0,1]
    # else :
    #     print 'No scale factor available for the current matching time:',mth_time
    #     sys.exit(-1)
    # find the scaling factor for time tau by comparing string
    tau_series = ['%g'%i for i in sfactor_list[:,0]] # convert to string
    try:
        tau_idx = tau_series.index('%g'%mth_time) #index (row number) of tau
    except ValueError:
        print 'getSfactor error!'
        print 'tau=%g'%mth_time+' is not in sfactor_list.dat!'
        sys.exit(-1)
    sfactor = sfactor_list[tau_idx,1] 
    return sfactor


def getEventAngle(angle_order):
    """ 
    Get the event angle for profile rotation. Read the last line of the file.
    """
    #load data file
    epxInitial_file_path = path.join(lmDirectory, 'data/Epx_initial.dat')
    epxInitial_data = np.loadtxt(epxInitial_file_path)
    epxInitial_data = epxInitial_data.reshape(-1)  #rewrite data into 1D array
    epxAngle_now = epxInitial_data[(2*angle_order-19)]  #python index ends with -1
                               # order 2 access the [-1, 5] or [-1, -2]
    return epxAngle_now


def runHydro(mth_time, hydro_dir, inital_phi):
    """
    run hydro code to get the freeze-out surface
    """
    #scale the profile
    norm_factor_now = "%.6f" %norm_factor
    tau0 = "%.2f" %mth_time
    angle_now = "%.6f" %inital_phi
    #run hydro
    print 'Hydro is running at tau0:', mth_time
    hydro_cmd = hydro_dir + "/./VISHNew.e >runlog.dat"  + " t0="+ tau0 \
                          + " factor=" + norm_factor_now  \
                          + " event_angle=" + angle_now                                                 
    retcode = call(hydro_cmd, shell=True, cwd=hydro_dir)
    # check if hydro excuted successfully
    if retcode ==0:
        print 'Hydro evolution completes!'
    else :
        print 'Hydro stops unexpectly!\n'
        sys.exit(-1)
    # check if hydro finishes
    runlog_file = path.join(hydro_dir,'runlog.dat')
    logEndLine = popen("tail -n 1 %s" %runlog_file).read() #
    logWord = logEndLine.split()[0]
    if logWord == 'Finished':
        success = True
    else:
        print 'Hydro did not run to the end!'
        success = False
    # backup the hydro run log
    shutil.copy(path.join(hydroDirectory, 'runlog.dat'),  \
        path.join(hydroResultDirectory, 'runlog.dat')) 
    return success



def cleanUpFolder(targetDir):
    """ 
    Delete all data files in the given directory. 
    """
    if path.exists(targetDir):
        try:
            call("rm -rf *", shell=True, cwd=targetDir)
        except OSError:
            pass # very likely the the folder is already empty
    else:
        makedirs(targetDir)



def runcodeShell():
    """
    run LandauMatching, Hydro, iS consequently and move the data files.
    """
    ##  move the inital KLN profiles to fs folder
    MCbackupDirectory = path.join(backupDir, 'MCevents')
    moveMCdata2Lm(MCbackupDirectory, lmInitDirectory)
    
    #   loop over the matching time
    for event_number in event_num_list:
        #   run LandauMatching and move the data files of one event at a specific 
        #   matching time to hydro inital conditions folder
        event_folder_tag = 'event_' + '%g' %event_number 
        runlm(event_number, sys_start_time, min(matching_time_list), \
            max(matching_time_list), matching_time_list[1]-matching_time_list[0])  

        for matching_time in matching_time_list:
            lmDataDirectory = path.join(lmDirectory, 'data/result/event_'
                                        + str(event_number)+"/"+"%g" %matching_time)
            cleanUpFolder(hydroInitialDirectory)
            moveLmData2Hydro(lmDataDirectory, hydroInitialDirectory, "%g" %matching_time)

            # run Hyrdo and move the results to iS
            cleanUpFolder(hydroResultDirectory)
            global norm_factor
            norm_factor = getSfactor(matching_time)
            angle_phi2 = getEventAngle(2)
            successCode = runHydro(matching_time, hydroDirectory, angle_phi2)
            if successCode == False:
            	break
            else:
            	print 'runcode: hydro does not finish successfully!'
            	print 'at event %d'%event_number + ', matching time %6.2f'%matching_time
            	sys.exit(-1)
            # backup the avg_points file
            backupFile = path.join(hydroResultDirectory, 'avg_points.dat')
            data_backup_dir = path.join(backupDir, event_folder_tag, 
                '%g'%matching_time, 'avg_points_corrected.dat')
            shutil.copy(backupFile,data_backup_dir)


if __name__ == "__main__":
    runcodeShell()

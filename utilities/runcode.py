#!/usr/bin/env python 

#   file: runcode.py
#
#   Author:   Jia Liu    <liu.2053@osu.edu>
#
#   This python script rename the initial profiles and move 
#   them to VISHNEW folder, then move the results folder of VISHNEW to iS
#
#   Revision history:
#   Apr.24,2014  Add function preProcessDecdat2File() to pre-process decdat2.dat file.
#               Now it calls doubleProjectPimn() to apply double projection to hydro dumped
#               pimn to make it traceless and orthogonal to u^mu.
#   Jan.31,2014  Match to total energy given by thermalization surface: decdat2.dat.
#               New 'Epx_initial.dat' format, getEventAngle() changes accordingly.
#   Jan.27.2014  Calculate and output dEdydphip for thermalization surface part and hydro
#                part.
#   Jan.24,2014  Adopted from the old version, converted to output dEdydphip table.
#   Sep.26.2013  Add more arguments to runlm, to control matching time range.
#   Sep.18.2013  Import angle to hydro running. Backup "Epx_inital.dat" and each hydro
#                "runlog.dat".
#   Aug.09.2013  Backup the hydro runlog.dat. Use it to check pi^munu regulation later.
#   Aug.06.2013  Modified to cluster usable version: including runsuperMC, storing 
#                data in a separate folder.
#   Aug.01.2013  Backup the data for each run in separate folders, inside
#                iS/backup/event_#/#.#. Then use matlab to extract data and analysis
#   May.30.2013  Incorporate scale factor for sqrt(s)=2760GeV, MC-KLN.
#   May.23.2013  For multiplicity's sake, divide the sFactor of hydro by 8  
#                to get the right thermal pion number: 200~220.
import sys
from os import path, getcwd, listdir, rename, makedirs, remove, stat, popen
import shutil #contains command to move files
from subprocess import call
import numpy as np
import dEcounters

event_num_list = range(1, 41)
events_total = len(event_num_list)
matching_time_list = [2,3,4,6,7,8,10]
pre_process_decdat2_file = True
sys_start_time = 0.01


angle_phi2 = 0.0
norm_factor = 1.0
sfactor_list = np.loadtxt('sfactor_log.dat')  #load in scale factor

rootDir = getcwd()  #define the root folder

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

#define the folder of iS
iSDirectory = path.join(rootDir, 'iS')
iSDataDirectory = path.join(iSDirectory, 'results')

backupDirName = rootDir.split('/')[-1] # get the name of the current node
backupDir = path.join('..','..', 'dataBase', backupDirName)




def runsuperMC(events_total):
    """
    Run superMC to get the inital KLN profile
    """
    superMC_cmd = superMCDirectory + '/./superMC.e '  \
        +' operation=1' + ' nev='+str(events_total)
    superMC_retcode = call(superMC_cmd, shell=True, cwd=superMCDirectory)
    if superMC_retcode ==0:
      print 'SuperMC completes!'
    else :
      print 'SuperMC stops unexpectly!'
      sys.exit(-1)    
    #store the parameters.dat for current run
    shutil.copy(path.join(superMCDirectory, 'parameters.dat'), path.join(superMCResultDirectory, 'parameters.dat'))


def backupMCdata(source_dir, dest_dir):
    """
    backup inital superMC generated data files
    """
    cleanUpFolder(dest_dir)
    for init_files in listdir(source_dir):
        src_file = path.join(source_dir, init_files)
        dst_file = path.join(dest_dir, init_files)
        shutil.copy(src_file, dst_file)
    print 'superMC profiles have been backed up!'  


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



def getInitialEcc2(file_path):
    """ Read in the initial eccentricity for one single event
    """
    ecc2_raw = np.loadtxt(file_path)
    global ecc2_initial
    ecc2_initial = ecc2_raw[-1,:] #just read the last line for the current event


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
    #find if the matching time is in the table
    temp = '%.8f'%mth_time
    mth_time = float(temp)
    if (sfactor_list[:,0]==mth_time).any()==True:
        return (sfactor_list[sfactor_list[:,0]==mth_time])[0,1]
    else :
        print 'No scale factor available for the current matching time:',mth_time
        sys.exit(-1) 

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
        pass
    else:
        print 'Hydro did not run to the end!'
        sys.exit(-1)
    # backup the hydro run log
    shutil.copy(path.join(hydroDirectory, 'runlog.dat'),  \
        path.join(hydroResultDirectory, 'runlog.dat')) 



def moveHydroData2iS(source_dir, dest_dir):
    """
    move the entire /results folder of hydro to iS
    """
    # move the hydro results to iS folder
    for data_files in listdir(source_dir):
      if data_files.startswith('.nfs'):  #ignore junk files
        continue
      src_file = path.join(source_dir, data_files)
      dst_file = path.join(dest_dir, data_files)
      shutil.move(src_file, dst_file)
    # print 'Hyro results have been put into iS folder!\n'   


def doubleProjectPimn(decdat2_line):
    """
    Applying double projection operator on hydro dumped pi^munu.
    Return the line with replaced pi^munu
    """
    gmn = np.mat([[1., 0., 0., 0.], \
        [0., -1, 0., 0.],   \
        [0., 0., -1., 0.],  \
        [0., 0., 0., -1.]]) # metric
    # get umu rank-1 matrix
    vx_fo = decdat2_line[4]
    vy_fo = decdat2_line[5]
    gamma_fo = 1.0/np.sqrt(1.0 - vx_fo**2- vy_fo**2 )
    ux_fo = decdat2_line[4]*gamma_fo
    uy_fo = decdat2_line[5]*gamma_fo    
    uz_fo = 0
    umu = np.mat([gamma_fo, ux_fo, uy_fo, uz_fo])
    # get pi^munu matrix
    pimn = np.zeros([1, 7])
    for i in range(1,7):
        pimn[0,i-1] = decdat2_line[i+12] # pi00, pi01, pi02, pi11, pi12, pi22
    pimn[0,-1] = decdat2_line[12] # pi33
    pimn_input = np.mat([[pimn[0, 0], pimn[0, 1], pimn[0, 2], 0.],    \
    [pimn[0, 1], pimn[0, 3], pimn[0, 4], 0.], \
    [pimn[0, 2], pimn[0, 4], pimn[0, 5], 0], \
    [0., 0., 0., pimn[0, 6]] ])

    # apply double projection operator on pi^munu
    pimn_ap = (gmn*gmn - umu.T*(umu*gmn))*(((gmn*gmn - umu.T*(umu*gmn))*pimn_input).T) \
        -1.0/3.0*(gmn - umu.T*umu)*((np.trace(pimn_input*gmn)-(umu*gmn)*pimn_input*(gmn*umu.T)))[0,0]
    # update pi^munu the decdat2_oneLine
    decdat2_line[12] = pimn_ap[3,3]
    decdat2_line[13] = pimn_ap[0,0]
    decdat2_line[14] = pimn_ap[0,1]
    decdat2_line[15] = pimn_ap[0,2]
    decdat2_line[16] = pimn_ap[1,1]
    decdat2_line[17] = pimn_ap[1,2]
    decdat2_line[18] = pimn_ap[2,2]
    return decdat2_line


def preProcessDecdat2File(source_folder):
    """
    Pre-process input decdat2 file, replacing the
    original decdat2.dat. Backup the original one.
    """
    decdat2_fileName = 'decdat2.dat' 
    decdat2_file = path.join(source_folder, decdat2_fileName)
    try:
        decdat2_input = np.loadtxt(decdat2_file)
    except:
        print 'preProcessPimn: cannot read file: '+decdat2_file
        sys.exit(-1)

    # process file
    if decdat2_input.ndim == 1: # only one line in decdat2.dat
        decdat2_oneLine = decdat2_input
        decdat2_output = doubleProjectPimn(decdat2_oneLine)
    else:
        decdat2_output = np.zeros((1, 20)) # prepare a zero line to accept the following data
        decdat2_file_length = decdat2_input.shape[0]
        for i in range(0, decdat2_file_length):
            decdat2_oneLine = decdat2_input[i,:]
            decdat2_oneLine_new = doubleProjectPimn(decdat2_oneLine)
            decdat2_output = np.concatenate((decdat2_output,np.mat(decdat2_oneLine_new)))
        decdat2_output = np.delete(decdat2_output, 0, axis=0) # delete the zero data in line 1

    # back up the old decdat2.dat and save the new one
    decdat2_backupFileName = 'decdat2_backup.dat'
    decdat2_backupFile = path.join(source_folder, decdat2_backupFileName)
    shutil.move(decdat2_file, decdat2_backupFile)
    if decdat2_input.ndim == 1:
        np.savetxt(decdat2_file, decdat2_output[None], fmt='%14.6E')
    else:
        np.savetxt(decdat2_file, decdat2_output, fmt='%14.6E')



def runiS(iS_dir):
    """
    run iS code to get the anisotropy flow
    """
    print 'Start to run iS...'
    #run iS to get vn of thermal particles
    iS_cmd = iS_dir + "/./iS.e >runlog.dat"
    retcode = call(iS_cmd, shell=True, cwd=iS_dir)
    if retcode ==0:
      print 'iS completes!'
    else :
      print 'iS stops unexpectly!\n'
      sys.exit(-1)



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



def v2_ecc2_dataWrite(e_num, mth_time, is_data_dir):
    """
    Record inital epsilon_2 and final v_2 data for each event and 
    every matching time.
    """
    #open data file and write
    ratio_file = open('v2_over_ecc2_test.dat', 'a')
#    ratio_file.write("#Event number          Matching Time          v2/ecc2  \n")
    ecc2_sevent_initial = ecc2_initial[2]
    #open inital eccentricity data file of hydro and read in ecc2
    ecc2_file_path = is_data_dir + "/ecc-init-r_power-2.dat"
    ecc2_data = np.loadtxt(ecc2_file_path)
    ecc2_now = ecc2_data[1,2]  #python index start with 0
                               #1 means order = 2; 2 means accessing length of ecc2
    #open final thermal pion v2inte table and read in v2
    v2_file_path = is_data_dir + "/thermal_211_integrated_vndata.dat"
    v2_data = np.loadtxt(v2_file_path)
    v2_now = v2_data[2,-1]   # 2 means order = 2; -1 means accessing the length of v2
    v2_ecc2_ratio = v2_now/ecc2_sevent_initial
    #write them in the data table
    ratio_file.write("%g %10.4f %16.8f %16.8f %16.8f %16.8f\n" \
               %(e_num, mth_time, ecc2_sevent_initial, ecc2_now, \
                v2_now, v2_ecc2_ratio))  
    ratio_file.close()


def backupEcc2v2Files(source_dir, dest_dir, tag):
    """ 
    Backup "ecc-init-r_power-2.dat", "thermal_211_integrated_vndata.dat"
    and "VISH2p1_tec.dat" for each run.
    """
    #identify the data files
    src_file = source_dir+'/ecc-init-r_power-2.dat'
    dst_file = dest_dir+'/ecc-init-r_power-2'+ tag + '.dat'
    #rename the data files
    src_file_new_name = source_dir+'/ecc-init-r_power-2'+ tag + '.dat'
    rename(src_file, src_file_new_name)
    if path.exists(dst_file):
       remove(dst_file)
    shutil.move(src_file_new_name, dst_file)

    #identify the data files
    src_file = source_dir+'/thermal_211_integrated_vndata.dat'
    dst_file = dest_dir+'/thermal_211_integrated_vndata'+ tag + '.dat'
    #rename the data files
    src_file_new_name = source_dir+'/thermal_211_integrated_vndata'+ tag + '.dat'
    rename(src_file, src_file_new_name)
    if path.exists(dst_file):
       remove(dst_file)
    shutil.move(src_file_new_name, dst_file)
    
    #identify the data files
    src_file = source_dir+'/VISH2p1_tec.dat'
    dst_file = dest_dir+'/VISH2p1_tec'+ tag + '.dat'
    #rename the data files
    src_file_new_name = source_dir+'/VISH2p1_tec'+ tag + '.dat'
    rename(src_file, src_file_new_name)
    if path.exists(dst_file):
       remove(dst_file)    
    shutil.move(src_file_new_name, dst_file)    


def backupEachRun(source_dir, dest_dir):
    """ 
    Backup the whole 'results' folder in 'iS' for each run.
    """
    cleanUpFolder(dest_dir)
    #
    src_files = listdir(source_dir)
    for file_name in src_files:
        full_file_name = path.join(source_dir, file_name)
        if (path.isfile(full_file_name)):
            shutil.copy2(full_file_name, dest_dir)
    print 'Backup Completes!'


def moveTempRun(source_dir, dest_dir):
    """ 
    move the temporary results
    """
    #
    src_files = listdir(source_dir)
    for file_name in src_files:
        src_name = path.join(source_dir, file_name)
        dst_name = path.join(dest_dir, file_name)

        if (path.isfile(src_name)):
            shutil.copy2(src_name, dst_name)
        elif (path.isdir(src_name)):
            if path.exists(dst_name):
                shutil.rmtree(dst_name)
            shutil.copytree(src_name, dst_name)
#    print 'Temporary backup Completes!'

def getTotaldEdy(Sigma_th_file, edInit_file, normFactor, Sigma_th_outputFileName, \
        Sigma_fo_inputFolder, Sigma_fo_outputFileName, targetFolder):
    """
    run two functions from dEcounters to get the total dEdy, also generate 
    dEdydphip tables for two parts. 
    """
    dEdy_therm = dEcounters.dEdydphipSigmaThermal(Sigma_th_file,edInit_file,\
        normFactor, targetFolder, Sigma_th_outputFileName)
    dEdy_fo = dEcounters.dEdydphipSigmaFO(Sigma_fo_inputFolder, targetFolder,\
        targetFolder, Sigma_fo_outputFileName)
    # dEdy_{fo}_{pion}/dEdy_{fo}_{total} = alpha_pion
    totaldEdy = dEdy_therm + dEdy_fo
    return totaldEdy  



def runcodeShell():
    """
    run LandauMatching, Hydro, iS consequently and move the data files.
    """
    ##  run superMC and move the inital KLN profiles to fs folder
    #runsuperMC(events_total)
    #MCbackupDirectory = path.join(backupDir, 'MCevents')
    #backupMCdata(superMCResultDirectory, MCbackupDirectory)
    #moveMCdata2Lm(superMCResultDirectory, lmInitDirectory)
    
    #   loop over the matching time
    for event_number in event_num_list:
        #   run LandauMatching and move the data files of one event at a specific 
        #   matching time to hydro inital conditions folder
        event_folder_tag = 'event_' + '%g' %event_number 
        runlm(event_number, sys_start_time, min(matching_time_list), \
            max(matching_time_list), matching_time_list[1]-matching_time_list[0])  

        # two-fold way     
        #runlm(event_number, sys_start_time, min(matching_time_list1), \
        #         max(matching_time_list1), matching_time_list1[1]-matching_time_list1[0])

        # backup the 1st run results
        # temp_directory = path.join(lmDirectory, 'data/temp')
        #currentRun_directory = path.join(lmDirectory, 'data/result/event_'+str(event_number))
        #cleanUpFolder(temp_directory)
        #moveTempRun(currentRun_directory, temp_directory)

        # combine the results of two lm runs
        #moveTempRun(temp_directory, currentRun_directory)
        #cleanUpFolder(temp_directory)

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
            runHydro(matching_time, hydroDirectory, angle_phi2)
            cleanUpFolder(iSDataDirectory)
            moveHydroData2iS(hydroResultDirectory, iSDataDirectory)

            #run iS to get the anisotrpy flow
            decdat_file = path.join(iSDataDirectory, 'decdat2.dat')
            if(stat(decdat_file).st_size==0) :
                print 'decdat.dat is empty!\n'
                dEdy_th_inputFile = path.join(hydroInitialDirectory, 'dEd2rdphip_kln.dat')
                ed_th_file = path.join(hydroInitialDirectory, 'ed_profile_kln.dat')
                dEdy_th_outputFileName = 'dEdydphipThermal.dat'
                dEcounters.dEdydphipSigmaThermal(dEdy_th_inputFile, ed_th_file,\
                    norm_factor, iSDataDirectory, dEdy_th_outputFileName)
            else:
                if pre_process_decdat2_file is True:
                    preProcessDecdat2File(iSDataDirectory)
                runiS(iSDirectory)
                # call functions in module dEcounters to get dEdydphip tables
                dEdy_th_inputFile = path.join(hydroInitialDirectory, 'dEd2rdphip_kln.dat')
                ed_th_file = path.join(hydroInitialDirectory, 'ed_profile_kln.dat')
                dEdy_th_outputFileName = 'dEdydphipThermal.dat'
                dEdy_fo_outputFileName = 'dEdydphipFO.dat'
                totaldEdyTest = getTotaldEdy(dEdy_th_inputFile, ed_th_file, norm_factor,\
                    dEdy_th_outputFileName, \
                    iSDirectory, dEdy_fo_outputFileName, iSDataDirectory)   

            #back up the data files for this run
            sub_folder_tag = '%g' %matching_time
            data_backup_dir = path.join(backupDir, event_folder_tag, sub_folder_tag)
            backupEachRun(iSDataDirectory, data_backup_dir)
            print 'Current run completed for event ' + str(event_number) \
            + ', matching time '+ str(matching_time) +'fm/c!\n'
    # backup initial eccentricity file
    Epx_initial_file = path.join(lmDirectory, 'data', 'Epx_initial.dat')
    data_backup_dir = path.join(backupDir)
    shutil.copy(Epx_initial_file, path.join(data_backup_dir, 'Epx_initial.dat')) 


if __name__ == "__main__":
    runcodeShell()

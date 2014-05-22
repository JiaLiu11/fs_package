#!/usr/bin/env python

# Jan 31, 2014  Combined ratio method and bi-section method.
# Find the rescaling factor of the initial energy-momentum tensor
#and store the results corresponding to correct factor. 
# Using Bi-Section method to find sfactor.
# 

from os import path,stat
from numpy import *
import runcode
import time
import shutil
import dEcounters

totaldEdyExpect = 1599.4  #from single-shot viscous hydro starting at 0.6fm/c, thermal pion+ num 210.
event_number = 99
mtime_list = linspace(1, 4, 4) 
pre_process_decdat2_file = True
#mtime_list = [10]
sfactorL = 0.01
sfactorR = 20.0

def getParamsFromCML(parserCML):
	""" get the matching time from command line
	"""


def fetchFinalMultiplicity(targetDir, targetFile):
	""" get the final multiplicity from VISHNEW results
	"""
	thermal_pion_file = path.join(targetDir, targetFile)
	if path.isfile(thermal_pion_file) is False:
		print 'No such file: ', thermal_pion_file
		exit(-1)
	else:
            v2data=loadtxt(thermal_pion_file)
            return v2data[0,1]

def getRescaleFactor(total_dEdy_now):
	""" find the final multiplicity given by hydro to the required multiplicity
	"""
	#Get the rescale factor
	rescale_factor = totaldEdyExpect/total_dEdy_now
	return rescale_factor

def getTotaldEdyOnly(dEdyd2rdphipFile, edFile, sfactor, dEdydphipthermFolder, \
		dEdydphipFileName, iSFolder, dEdydphipHydroFolder, dEdydphipHydroFileName):
    """
    run two functions from dEcounters to get the total dEdy.
    """
    dEdy_therm = dEcounters.dEdydphipSigmaThermal(dEdyd2rdphipFile, edFile, sfactor, \
		dEdydphipthermFolder, dEdydphipFileName)
    iSDataFolder = path.join(iSFolder, 'results')
    dEdy_fo = dEcounters.dEdydphipSigmaFO(iSFolder, iSDataFolder, dEdydphipHydroFolder, dEdydphipHydroFileName)

    totaldEdy = dEdy_therm + dEdy_fo
    return totaldEdy


def adjustSfactorShell():
	""" use as main program
	"""
	sfactor_log = open('sfactor_log.dat', 'a')
	sfactor_log.write("#Matching Time         sfactor                 totaldEdy  \n")

        for matching_time in mtime_list:
    	    #Prepare the data for hydro run
    	    rescale_factor = 6.0
    	    rescale_factor_used = 6.0
            lmDataDirectory = path.join(runcode.lmDirectory, 'data/result/event_' \
                                        + str(event_number)+"/"+ "%g" %matching_time)
            runcode.cleanUpFolder(runcode.hydroInitialDirectory)
            runcode.moveLmData2Hydro(lmDataDirectory, runcode.hydroInitialDirectory, "%g" %matching_time)
	    xl = sfactorL
	    xr = sfactorR            
	    while 1:
		# use bisect at late matching
		if(matching_time>=7.9):
			runcode.norm_factor = (xl+xr)/2.0
		else:
			runcode.norm_factor = rescale_factor

		runcode.runHydro(matching_time, runcode.hydroDirectory, 0.0)
		runcode.cleanUpFolder(runcode.iSDataDirectory)
		runcode.moveHydroData2iS(runcode.hydroResultDirectory, runcode.iSDataDirectory)

		#run iS to get the total energy
		decdat_file = path.join(runcode.iSDataDirectory, 'decdat2.dat')
		if(stat(decdat_file).st_size==0) :  
			print 'decdat.dat is empty! Only calculate thermalization surface. \n'	
			dEdphipFile = path.join(runcode.hydroInitialDirectory, 'dEdydphip_kln.dat')		
			totaldEdyTest =dEcounters.dEdyTotalInital(dEdphipFile, \
				runcode.norm_factor)		
		else:
			if pre_process_decdat2_file is True:
				runcode.preProcessDecdat2File(runcode.iSDataDirectory)
				print 'adjustsfactor: use preprocessed decdat2 file!'
	   		runcode.runiS(runcode.iSDirectory)
			#get the current total energy
			dEdyd2rdphipFile = path.join(runcode.hydroInitialDirectory, 'dEd2rdphip_kln.dat')
			edFile = path.join(runcode.hydroInitialDirectory, 'ed_profile_kln.dat')

			totaldEdyTest = getTotaldEdyOnly(dEdyd2rdphipFile, edFile, runcode.norm_factor, \
				runcode.iSDataDirectory, 'dEdydphip_therm.dat',\
				runcode.iSDirectory, runcode.iSDataDirectory, 'dEdydphip_fo.dat')   

	  	#Keep a log	
		sfactor_log.write("%6.2f       %20.8f        %20.4f\n"   \
	   		       %(matching_time, runcode.norm_factor, totaldEdyTest))
	   	sfactor_log.flush()

        	#conditions for ending the loop
		if( abs(totaldEdyTest-totaldEdyExpect) < 10):#debug
		    sfactor_log.write("#Finally the factor:\n")
	   	    sfactor_log.write("%6.2f       %20.8f        %20.4f\n"   \
	                %(matching_time, runcode.norm_factor, totaldEdyTest))
	   	    #backup the run files
	   	    event_folder_tag = 'event_' + '%g' %event_number 
                    sub_folder_tag = '%g' %matching_time
		    data_backup_dir = path.join(runcode.backupDir, \
          	       event_folder_tag, sub_folder_tag)
		    runcode.backupEachRun(runcode.iSDataDirectory, data_backup_dir)
		    print 'Current run completed for event ' + str(event_number) 
   #          runlog_file = path.join(runcode.iSDataDirectory, 'runlog.dat')
			# runlog_newName = 'runlog_' + str(matching_time) +'_'    \
			#         + '%g' %rescale_factor_used + '.dat'
			# runlog_newFile = path.join(runcode.rootDir, runlog_newName)
			# shutil.move(runlog_file, runlog_newFile)
		    break	

	   	#get the used rescale factor for the next run
	   	rescale_factor_used = runcode.norm_factor
	   	if(matching_time>=7.9):
			if(totaldEdyTest-totaldEdyExpect<0):
				xl = rescale_factor_used
				xr = xr
			elif(totaldEdyTest-totaldEdyExpect>0):
				xl = xl
				xr = rescale_factor_used
		else:   	
			rescale_factor = getRescaleFactor(totaldEdyTest) * rescale_factor_used
			                	
	sfactor_log.close()

if __name__ == "__main__":
    adjustSfactorShell()

#!/usr/bin/env python

# Author: Jia Liu

# Purpose: 
#	Run fs+VISHNEW+iS+resDecay for different combination
#	of tau_s, eta/s, and T_dec.
#
# Implementation:
#	Specify the switching time for one node, then run two 
#	parameter search in adjust sfactor mode. After finding 
#	the correct scaling factor, run resonance decay 
#	and save the files to database. 
#

# libraries
from os import path,stat, popen
from subprocess import call
import numpy as np
import runcode, dEcounters
import sys, shutil

# parameters
taus_list = np.linspace(0.4,2.0,9) 
eta_s_list= np.linspace(0.08, 0.4, 9)
tdec_list = np.linspace(100, 160, 7) #MeV
edec_list = np.array([0.097561,0.135765,0.1822,0.237036,
	0.300612,0.373435,0.456113]) # GeV/fm^3, EOS 95pceV0

# folder structure
rootDir = path.abspath("..") #assume this file is placed under utilities/

# run parameters
node_index = rootDir.split('/')[-1]  # get the name of current node
taus_index = int(node_index.split('e')[-1])-1 # index of tau_s for current folder
matching_time= taus_list[taus_index]
backupDir_currentNode = path.join('..','..','..', 'dataBase', 'taus_%g'%matching_time)

# parameters for scaling factor search
totaldEdyExpect = 1599.4  #from single-shot viscous hydro starting at 0.6fm/c, thermal pion+ num 210.
event_number = 99
pre_process_decdat2_file = True
sfactorL = 0.01
sfactorR = 20.0
rescale_factor_guess = 9.0 #initial guess of rescaling factor


def runlm(event_num, tau0, taumin, taumax, dtau, lmDirectory):
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


def runHydro(mth_time, norm_factor, hydro_dir, inital_phi, etas, edec):
    """
    run hydro code to get the freeze-out surface
    """
    #scale the profile
    norm_factor_now = "%.6f" %norm_factor
    tau0 = "%.2f" %mth_time
    angle_now = "%.6f" %inital_phi
    etas_now = "%.g" %etas
    edec_now = "%.g" %edec
    #run hydro
    print 'Hydro is running at tau0:', mth_time
    hydro_cmd = hydro_dir + "/./VISHNew.e >runlog.dat"  + " t0="+ tau0 \
                          + " factor=" + norm_factor_now  \
                          + " event_angle=" + angle_now   \
                          + " edec=" + edec_now \
                          + " viscousc="+  etas_now                                           
    retcode = call(hydro_cmd, shell=True, cwd=hydro_dir)
    # check if hydro excuted successfully
    if retcode ==0:
        pass
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
    shutil.copy(path.join(runcode.hydroDirectory, 'runlog.dat'),  \
        path.join(runcode.hydroResultDirectory, 'runlog.dat')) 
    return success

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

def parameterSearchShell():
	# get free-streamed data
	runlm(99, 0.01, matching_time, matching_time+0.01, 0.01, runcode.lmDirectory)

	# move lm data to hydro
	lmDataDirectory = path.join(runcode.lmDirectory, 'data/result/event_' \
	                            + str(event_number)+"/"+ "%g" %matching_time)
	runcode.cleanUpFolder(runcode.hydroInitialDirectory)
	runcode.moveLmData2Hydro(lmDataDirectory, runcode.hydroInitialDirectory, \
		"%g" %matching_time)

	sfactor_log = open(path.join(rootDir,'sfactor_log.dat'), 'a+')
	sfactor_log.write("##Matching Time             eta/s                 Tdec                  sfactor                 totaldEdy \n")	

	# start to loop over two parameter space
	rescale_factor_used = rescale_factor_guess
	for eta_s in eta_s_list:
		for tdec_idx in range(len(tdec_list)):
			tdec = tdec_list[tdec_idx]
			edec = edec_list[tdec_idx]
			rescale_factor = rescale_factor_used
			xl = sfactorL
			xr = sfactorR 
			print 'Running at taus=%g, eta/s=%.3f, Tdec=%.3f' \
				%(matching_time, eta_s, tdec)           
			while 1:
				# use bisect at late matching
				if(matching_time>=7.9):
					norm_factor = (xl+xr)/2.0
				else:
					norm_factor = rescale_factor

				runHydro(matching_time, norm_factor, runcode.hydroDirectory, 0.0, 
					eta_s, edec) 
				runcode.cleanUpFolder(runcode.iSDataDirectory)
				runcode.moveHydroData2iS(runcode.hydroResultDirectory, runcode.iSDataDirectory)

				#run iS to get the total energy
				decdat_file = path.join(runcode.iSDataDirectory, 'decdat2.dat')
				if(stat(decdat_file).st_size==0) :  
					print 'decdat.dat is empty! Only calculate thermalization surface. \n'	
					dEdphipFile = path.join(runcode.hydroInitialDirectory, 'dEdydphip_kln.dat')		
					totaldEdyTest =dEcounters.dEdyTotalInital(dEdphipFile, \
						norm_factor)		
				else:
					if pre_process_decdat2_file is True:
						runcode.preProcessDecdat2File(runcode.iSDataDirectory)
						print 'adjustsfactor: use preprocessed decdat2 file!'
					runcode.runiS(runcode.iSDirectory)
					#get the current total energy
					dEdyd2rdphipFile = path.join(runcode.hydroInitialDirectory, 'dEd2rdphip_kln.dat')
					edFile = path.join(runcode.hydroInitialDirectory, 'ed_profile_kln.dat')

					totaldEdyTest = getTotaldEdyOnly(dEdyd2rdphipFile, edFile, norm_factor, \
						runcode.iSDataDirectory, 'dEdydphipThermal.dat',\
						runcode.iSDirectory, runcode.iSDataDirectory, 'dEdydphipFO.dat')   

				#conditions for ending the loop
				if( abs(totaldEdyTest-totaldEdyExpect) < 10):
					sfactor_log.write("%6.2f       %20.8f    %20.8f		%20.8f    %20.4f\n"   \
						%(matching_time, eta_s, tdec, norm_factor, totaldEdyTest))

					# run resonance decay
					runRes_cmd = runcode.iSDirectory+'/./resonance.e >reslog.dat'
					reso_code = call(runRes_cmd, shell=True, cwd=runcode.iSDirectory)
					if reso_code == 0:
					    pass
					else:
					    print 'parameterSearch: resonance decay failed!'
					    sys.exit(-1)	
					# run iInteSp
					inteSp_cmd = runcode.iSDirectory+'/./iInteSp.e'
					inteSp_code = call(inteSp_cmd, shell=True, cwd=runcode.iSDirectory)
					if inteSp_code == 0:
					    pass
					else:
					    print 'parameterSearch: iInteSp failed!'
					    sys.exit(-1)			
					#backup the run files
					backupDir = path.join(backupDir_currentNode, 'etas_%g'%eta_s, 
						'tdec_%g'%tdec)
					runcode.backupEachRun(runcode.iSDataDirectory, backupDir)
					print 'Current run completed for eta/s=%.3f, Tdec=%.3f' \
						%(eta_s, tdec)
					#get the used rescale factor for the matching time run
					rescale_factor_used = norm_factor
					break	
				# if current scaling facotor is not correct, get the guess for next run
				#Keep a log	
				sfactor_log.write("%6.2f       %20.8f    %20.8f		%20.8f    %20.4f\n"   \
					%(matching_time, eta_s, tdec, norm_factor, totaldEdyTest))
				sfactor_log.flush()
				if(matching_time>=7.9):
					if(totaldEdyTest-totaldEdyExpect<0):
						xl = norm_factor
						xr = xr
					elif(totaldEdyTest-totaldEdyExpect>0):
						xl = xl
						xr = norm_factor
				else:   	
					rescale_factor = getRescaleFactor(totaldEdyTest) * norm_factor

	sfactor_log.close()
	print "All parameters complete!"

if __name__ == "__main__":
    parameterSearchShell()


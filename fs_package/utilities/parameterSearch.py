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

# folder structure
rootDir = path.abspath("..") #assume this file is placed under utilities/
tableLocation = path.join(rootDir, "tables")

# run parameters
node_name = rootDir.split('/')[-1]  # get the name of current node
node_index = int(node_name.split('e')[-1]) # index of tau_s for current folder
backupDir_currentNode = path.join(rootDir, 'localdataBase', node_name)
number_of_nodes = 20 # total number of nodes
resumeMode = True
skipEvent = True

# parameters for scaling factor search
totaldEdyExpect = 1547  #from single-shot viscous hydro starting at 0.6fm/c, thermal pion+ num 210. 1547 for glb
event_number = 99 # 99: for smooth event which maximizes epsilon_2
pre_process_decdat2_file = True
sfactorL = 0.01
sfactorR = 20.0
rescale_factor_guess = 20 #initial guess of rescaling factor, 20 for glb

def findResumeRunNumber(param_log_fileName):
	"""
	find the number of events has been finished by counting how many lines in
	the parameter search log file of the current node.
	"""
	finished_events = 0
	if(path.isfile(param_log_fileName)):
		try:
			log_data = np.loadtxt(param_log_fileName)
		except:
			print "Empty parameter search log file!"
		if(log_data.ndim==1):
			finished_events = 0 if(log_data.shape[0]==0) else 1
		else:
			finished_events = log_data.shape[0] # a line represents a finished parameter search
			print "Resume run from event: %d"%(finished_events+1)
	else:
		print "Parameter log file does not exists!"
	return finished_events


def extractParameterList(infile_name, outfile_name, node_idx):
	"""
	Extract the parameter list for current node
	Return a numpy array
	"""
	params_list_data = np.loadtxt(infile_name)
	# subset data
	total_lines = params_list_data.shape[0]
	node_lines  = int(total_lines/number_of_nodes)
	params_selected = params_list_data[(node_idx-1)*node_lines:(node_idx*node_lines),:]
	np.savetxt(outfile_name, params_selected, fmt='%g',delimiter='\t')


def runlm(event_num, tau0, taumin, taumax, dtau, lmDirectory):
    """
    Run Landau matching code for a single event
    """
    lm_cmd = lmDirectory + "/./lm.e >runlog_lm.dat" + ' event_mode='+str(event_num) \
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
    tau0 = "%g" %mth_time
    angle_now = "%.6f" %inital_phi
    etas_now = "%g" %etas
    edec_now = "%g" %edec
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
  #   dEdy_therm = dEcounters.dEdydphipSigmaThermal(dEdyd2rdphipFile, edFile, sfactor, \
		# dEdydphipthermFolder, dEdydphipFileName) # for no pre-equilibrium
    iSDataFolder = path.join(iSFolder, 'results')
    dEdy_fo = dEcounters.dEdydphipSigmaFO(iSFolder, iSDataFolder, dEdydphipHydroFolder, dEdydphipHydroFileName)

    totaldEdy = dEdy_fo #dEdy_therm + dEdy_fo
    return totaldEdy

def storeParamSearchResults(event_id, data_folder, filename):
	"""
	Collect and write the parameter search result to log file.
	Result: write a new line contains: 
		event_id, taus_now, etas_now, tdec_now, v2ch_now, v3ch_now, meanpT_now
	"""
	# pT gaussian points and weights
	pT_new = np.array([0.33401871,0.37228685,0.44019231,0.53615897, \
	             0.65793937,0.80267954,0.96698706,1.14701076, \
	             1.33853107,1.53705894,1.73794106,1.93646893, \
	             2.12798924,2.30801294,2.47232046,2.61706063, \
	             2.73884103,2.83480769,2.90271315,2.94098129]);
	pT_weight = np.array([2.31183844e-02,5.32893766e-02,8.22570634e-02, \
	              1.09300723e-01,1.33783282e-01,1.55130323e-01, \
	              1.72841338e-01,1.86501143e-01,1.95789545e-01, \
	              2.00488821e-01,2.00488821e-01,1.95789545e-01, \
	              1.86501143e-01,1.72841338e-01,1.55130323e-01, \
	              1.33783282e-01,1.09300723e-01,8.22570634e-02, \
	              5.32893766e-02,2.31183844e-02]);
	# get parameter list
	params_now = np.loadtxt(path.join(data_folder, 'params.dat'))
	taus_run, etas_run, tdec_run, edec_run = params_now
	# get flow anisotropy
	vn_ch_data = np.loadtxt(path.join(data_folder, 'Charged_ptcut0510_eta_integrated_vndata.dat'))
	v2_ch = vn_ch_data[2, -1]
	v3_ch = vn_ch_data[3, -1]
	# get pion mean pt at pt range (0.3, 3) GeV/c
	pion_spectra_data = np.loadtxt(path.join(data_folder, 'pion_p_vndata.dat'))
	pt_array = pion_spectra_data[:, 0]
	dndpt_array = pion_spectra_data[:, 2] # dN/(2pi pT dpT)
	dndpt_interped = np.exp(np.interp(pT_new, pt_array, np.log(dndpt_array)))
	dn_interped = dndpt_interped*2.0*np.pi*pT_new
	meanpt_pion = sum(dn_interped*pT_new*pT_weight)/sum(dn_interped*pT_weight)
	# get proton mean pt at pt range (0.3, 3) GeV/c
	proton_spectra_data = np.loadtxt(path.join(data_folder, 'proton_vndata.dat'))
	pt_array = proton_spectra_data[:, 0]
	dndpt_array = proton_spectra_data[:, 2] # dN/(2pi pT dpT)
	dndpt_interped = np.exp(np.interp(pT_new, pt_array, np.log(dndpt_array)))
	dn_interped = dndpt_interped*2.0*np.pi*pT_new
	meanpt_proton = sum(dn_interped*pT_new*pT_weight)/sum(dn_interped*pT_weight)
	# write to file
	param_log = open(filename, 'a+')
	param_log.write("%d \t %12.6f \t %12.6f \t %12.6f \t %12.4e \t %12.4e \t %12.6f \t %12.6f\n"\
		%(event_id, taus_run, etas_run, tdec_run, v2_ch, v3_ch, meanpt_pion, meanpt_proton))
	param_log.close()


def parameterSearchShell():
	# get search parameters for current nodes
	params_list_source = path.join(tableLocation, "params_list_glb.dat") #"params_list_glb.dat" for glb
	params_list_current= path.join(rootDir, "params_list_"+
		node_name+".dat")
	extractParameterList(params_list_source, params_list_current, node_index)
	params_list_currentNode = np.loadtxt(params_list_current)

	rescale_factor_used = rescale_factor_guess # start value for sfactor search
	sfactor_log = open(path.join(rootDir,'sfactor_log.dat'), 'a+')
	sfactor_log.write("##Matching Time \t eta/s \t Tdec \t sfactor  \t totaldEdy \n")

	# store parameter search result
	paramSearchLogFile = open(path.join(rootDir, 'param_search_log.dat'),'a+') 
	paramSearchLogFile.write("#   tau_s \t eta/s \t tdec \t v2_ch  \t v3_ch \t meanpT_pion \t meanpT_proton \n")
	paramSearchLogFile.close()
	param_search_log = path.join(rootDir, 'param_search_log.dat')

	# store skipped parameters
	paramSearchSkipLogFile = open(path.join(rootDir, 'param_search_skipLog.dat'),'a+')
	paramSearchSkipLogFile.write("# event_num \t tau_s \t eta/s \t tdec \t e_dec \n")
	paramSearchSkipLogFile.close()
	param_search_skipLog = path.join(rootDir, 'param_search_skipLog.dat')

	# start event number
	startFromNum = 0
	if(resumeMode==True):
		startFromNum = findResumeRunNumber(param_search_log)

	for i in range(startFromNum, params_list_currentNode.shape[0]):
		matching_time, eta_s, tdec, edec = params_list_currentNode[i,:]
		# skip a parameter combination
		if(skipEvent==True and matching_time<0.1):
			print 'Skipping taus=%g, eta/s=%.3f, Tdec=%.3f ... ' \
				%(matching_time, eta_s, tdec)  
			# keep a log of what parameter skipped
			skipLog = open(param_search_skipLog, 'a+')
			skipLog.write("%d \t %g \t %g \t %g \t %g\n"\
				%(i+1, matching_time, eta_s, tdec, edec))
			skipLog.flush()
			skipLog.close()
			continue
		print 'Running at taus=%g, eta/s=%.3f, Tdec=%.3f ... ' \
			%(matching_time, eta_s, tdec)  
		# # get free-streamed data # for no pre-equilibrium
		# runlm(event_number, 0.01, matching_time, matching_time+0.01, 0.01, runcode.lmDirectory)
		# lmDataDirectory = path.join(runcode.lmDirectory, 'data/result/event_' \
		#                             + str(event_number)+"/"+ "%g" %matching_time)
		# runcode.cleanUpFolder(runcode.hydroInitialDirectory)
		# runcode.moveLmData2Hydro(lmDataDirectory, runcode.hydroInitialDirectory, \
		# 	"%g" %matching_time)
		# lmDataDirectory_upper = path.join(runcode.lmDirectory, 'data/result/event_' \
		# 	+ str(event_number))
		# runcode.cleanUpFolder(lmDataDirectory_upper)

		# start paramter search	
		rescale_factor = rescale_factor_used
		xl = sfactorL
		xr = sfactorR          
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
				# dEdphipFile = path.join(runcode.hydroInitialDirectory, 'dEdydphip_kln.dat') # for no prequilibrium		
				# totaldEdyTest =dEcounters.dEdyTotalInital(dEdphipFile, \
				# 	norm_factor)	
				print "empty decdat2.dat!"
				sys.exit(-1)	
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
				sfactor_log.write("%8.4f       %20.8f    %20.8f		%20.8f    %20.4f\n"   \
					%(matching_time, eta_s, tdec, norm_factor, totaldEdyTest))
				sfactor_log.flush()
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

				#backup the run files by zip
				backupDir = path.join(backupDir_currentNode, 'run_%d'%(i+1))
				runcode.backupEachRun(runcode.iSDataDirectory, backupDir)
				#write the run parameters to backup folder
				params_backup_file = open(path.join(backupDir, "params.dat"), "w")
				params_backup_file.write("# tau_s     eta_s     t_dec      e_dec\n")
				params_backup_file.write("%g 	%g 	%g 	%g\n" \
					%(matching_time, eta_s, tdec, edec))
				params_backup_file.close() 
				# write down parameter search result
				storeParamSearchResults(event_number, backupDir, param_search_log)
				# zip the backupdir and delete source files
				zip_cmd = "zip -r -q -j -m "+ "run_%d.zip "%(i+1)+ backupDir 
				zip_code = call(zip_cmd, shell=True, cwd=backupDir_currentNode)
				if zip_code == 0:
					pass
				else:
					print 'zip fails!'
					sys.exit(-1)
				print 'Current run completed (%d / %d ) for eta/s=%.3f, Tdec=%.3f\n' \
					%(i+1, params_list_currentNode.shape[0], eta_s, tdec)

				#get the used rescale factor for the matching time run
				rescale_factor_used = norm_factor
				break	
			# if current scaling facotor is not correct, get the guess for next run
			#Keep a log	
			sfactor_log.write("%8.4f       %20.8f    %20.8f		%20.8f    %20.4f\n"   \
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


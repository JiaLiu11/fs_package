#!/usr/bin/env python

# purpose: rerun iInteSp to recalculate the all charged particle flow
#		   anisotropy at |eta|<0.5 in accordance with ATLAS data.


# libraries
from os import path,stat, popen
from subprocess import call
import numpy as np
import runcode, dEcounters
import sys, shutil

# run parameters
event_num = 99 # for v2
totalEvents_num = 1300
nodes_num = 13
enable_skip = True

# folder structure
rootDir = path.abspath("..") #assume this file is placed under utilities/
node_name = rootDir.split('/')[-1]  # get the name of current node
node_index = int(node_name.split('e')[-1]) # index of events for current folder
source_folder = path.join(path.abspath("../../../"), 'zipped_1300_runs')
backupDir= path.join(rootDir, 'localdataBase', node_name)
runcode.cleanUpFolder(backupDir)

# calculate run range
runs_perNode = int(totalEvents_num/nodes_num)
runs_begin = (node_index-1)*runs_perNode+1
runs_end   = node_index*runs_perNode
currentRuns_range= range(runs_begin, runs_end+1)



def storeMultiplicityResults(event_id, data_folder, param_log):
	"""
	Collect and write the parameter search result to log file.
	Result: write a new line contains: 
		event_id, taus_now, etas_now, tdec_now, multiplicity
	"""
	# get parameter list
	params_now = np.loadtxt(path.join(data_folder, 'params.dat'))
	taus_run, etas_run, tdec_run, edec_run = params_now
	# get flow anisotropy
	dN_ch_rawdata = np.loadtxt(path.join(data_folder, 'Charged_eta_integrated_vndata.dat'))
	dn_ch = dN_ch_rawdata[0, 1] 
	# write to file
	param_log.write("%12.6f \t %12.6f \t %12.6f \t %12.6f \n"\
		%(taus_run, etas_run, tdec_run, dn_ch))
	param_log.flush()


def reCalMultiplicity():
	paramSearchLogFile = open(path.join(rootDir, 'param_search_multiplicity_log.dat'),'a+') # store parameter search result
	paramSearchLogFile.write("#    tau_s             eta/s         tdec           dNdeta\n")
	print "Begin to recalculate charged multiplicity...\n"
	print "current run from %d to %d \n"%(runs_begin, runs_end)

	for irun in currentRuns_range:
		# copy data to iS results
		runcode.cleanUpFolder(runcode.iSDataDirectory)
		run_zip_source = path.join(source_folder, 'run_%d.zip'%irun)
		if(path.isfile(run_zip_source)==False):
			if(enable_skip==True):
				print "*** skip run: run_%d.zip"%irun
				continue
			else:
				print "No such run: run_%d.zip"%irun
				sys.exit(-1)
		else:
			shutil.copy(run_zip_source, runcode.iSDataDirectory)

		# decompress zip 
		unzip_cmd = "unzip -q -j run_%d.zip && rm -f run_%d.zip "%(irun,irun)
		unzip_code = call(unzip_cmd, shell=True, cwd=runcode.iSDataDirectory)
		if unzip_code == 0:
			pass
		else:
			print 'zip fails!'
			sys.exit(-1)

		# run iInteSp
		inteSp_cmd = runcode.iSDirectory+'/./iInteSp.e'
		inteSp_code = call(inteSp_cmd, shell=True, cwd=runcode.iSDirectory)
		if inteSp_code == 0:
		    pass
		else:
		    print 'parameterSearch: iInteSp failed!'
		    sys.exit(-1)	

		# dump results
		storeMultiplicityResults(event_num, runcode.iSDataDirectory, paramSearchLogFile)
		# compress and store the results
		zip_cmd = "zip -q -j -m run_%d.zip ./* "%irun
		zip_code = call(zip_cmd, shell=True, cwd=runcode.iSDataDirectory)
		if zip_code == 0:
			pass
		else:
			print 'zip fails!'
			sys.exit(-1)
		zip_source = path.join(runcode.iSDataDirectory, 'run_%d.zip'%irun)
		shutil.copy(zip_source, backupDir)
		print 'run %d finished!'%irun
	paramSearchLogFile.close()
	print 'All jobs finished!'

if __name__ == "__main__":
    reCalMultiplicity()

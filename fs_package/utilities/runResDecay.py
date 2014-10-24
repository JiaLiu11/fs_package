#!/usr/bin/env python

#   file: runResDecay.py
#	purpose: 
#   Author:   Jia Liu    <liu.2053@osu.edu>
#   History:
#   Jul 28, 2014     First version.

import sys, shutil, glob
from os import path, stat, getcwd, stat, rename, remove, mkdir
import numpy as np
from subprocess import call

# run parameters
event_list = range(1,41)
#matchingTime_list = np.linspace(0.4,2.0,9)
matchingTime_list = [0.6, 0.8, 1, 1.2, 2]

# file locations and folders locations
rootDir = path.abspath('..')  #root directory for current node
database_location = path.join(path.abspath('../../..'),'dataBase')
is_location = path.join(rootDir, 'iS')
node_name = rootDir.split('/')[-1]

# files needed by resonance decay
dN_file_iSfolder = path.join(is_location,'results')
dN_file_target = path.join(is_location,'results','dN_ptdptdphidy.dat')
for event_num in event_list:
	for tau_s in matchingTime_list:
		# copy dNdyptdptdphi file from database to current iS if it exists
		dN_file_sourceFolder = path.join(database_location, node_name, 
			"event_%d"%event_num, "%g"%tau_s)
		dN_file_sourceFile = path.join(database_location, node_name, 
			"event_%d"%event_num, "%g"%tau_s, "dN_ptdptdphidy.dat")
		if(path.isfile(dN_file_sourceFile)):
			print "**** tau_s = %g is calculating"%tau_s
			shutil.rmtree(dN_file_iSfolder)# clear iS results folder
			mkdir(dN_file_iSfolder)
			shutil.copy(dN_file_sourceFile, dN_file_target)
			# run resonance decay
			reso_cmd = is_location+'/./resonance.e >log.dat'
			reso_code = call(reso_cmd, shell=True, cwd=is_location)
			if reso_code == 0:
			    pass
			else:
			    print 'runResDecay: resonance decay failed!'
			    sys.exit(-1)
			# run iInteSp
			inteSp_cmd = is_location+'/./iInteSp.e'
			inteSp_code = call(inteSp_cmd, shell=True, cwd=is_location)
			if inteSp_code == 0:
			    pass
			else:
			    print 'runResDecay: iInteSp failed!'
			    sys.exit(-1)
			# backup files from iS results folder to database
			for files in glob.glob(path.join(dN_file_iSfolder, '*.*')):
				shutil.copy(files, dN_file_sourceFolder)
		else:
			print "**** tau_s = %g is skipped"%tau_s
			continue
	print "runResDecay: event %d"%event_num+" processed!"
print "All %d events processed!"%len(event_list)

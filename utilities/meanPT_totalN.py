#!/usr/bin/env python

#   file: meanPT_totalN.py
#
#   Author:   Jia Liu    <liu.2053@osu.edu>
#   History:
#   Jul 16, 2014     First version.

import sys, shutil
from os import path, stat, getcwd, stat, rename, remove, mkdir
import numpy as np
from subprocess import call

# run parameters
event_list = range(1,41)
matchingTime_list = np.linspace(0.4,1.8,8)

# file locations
# folders locations
rootDir = path.abspath('..')  #root directory for current node
database_location = path.join(path.abspath('../../..'),'dataBase')
is_location = path.join(rootDir, 'iS')
node_name = rootDir.split('/')[-1]

# files needed by resonance decay
dN_file_iSfolder = path.join(is_location,'results')
dN_file_target = path.join(is_location,'results','dN_ptdptdphidy.dat')
for event_num in event_list:
	for tau_s in matchingTime_list:
		# copy dNdyptdptdphi file from database to current iS
		dN_file_source = path.join(database_location, node_name, 
			"event_%d"%event_num, "%g"%tau_s, "dN_ptdptdphidy.dat")
		shutil.rmtree(dN_file_iSfolder)# clear iS results folder
		mkdir(dN_file_iSfolder)
		shutil.copy(dN_file_source, dN_file_target)

		# run resonance decay
		reso_cmd = is_location+'/./resonance.e >log.dat'
		reso_code = call(reso_cmd, shell=True, cwd=is_location)
		if reso_code == 0:
		    pass
		else:
		    print 'meanPT_totalN: resonance decay failed!'
		    sys.exit(-1)
		# run iInteSp
		inteSp_cmd = is_location+'/./iInteSp.e'
		inteSp_code = call(inteSp_cmd, shell=True, cwd=is_location)
		if inteSp_code == 0:
		    pass
		else:
		    print 'meanPT_totalN: iInteSp failed!'
		    sys.exit(-1)
		# backup files from iS results folder to database
		for filename in glob.glob(path.join(dN_file_iSfolder, '*.*')):
    		shutil.copy(filename, dN_file_source)
	print "event %d"%event_num+" processed!"
print "All %d events processed!"%len(event_list)

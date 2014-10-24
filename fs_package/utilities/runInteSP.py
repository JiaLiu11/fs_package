#!/usr/bin/env python

#   file: runInteSP.py
#	purpose: 
#   Author:   Jia Liu    <liu.2053@osu.edu>
#   History:
#   Oct 06, 2014     First version.

import sys, shutil, glob
from os import path, stat, getcwd, stat, rename, remove, mkdir
import numpy as np
from subprocess import call

# run parameters
taus_list = np.linspace(0.4,2.0,9) 
eta_s_list= np.linspace(0.08, 0.4, 9)
tdec_list = np.linspace(100, 160, 7) #MeV

# file locations and folders locations
rootDir = path.abspath("..")
node_index = rootDir.split('/')[-1]  # get the name of current node
taus_index = int(node_index.split('e')[-1])-1 # index of tau_s for current folder
matching_time= taus_list[taus_index]
backupDir_currentNode = path.join('..','..','..', 'dataBase', 'taus_%g'%matching_time)
is_location = path.join(rootDir, 'iS')

# files needed by resonance decay
dN_file_iSfolder = path.join(is_location,'results')
dN_file_target = path.join(is_location,'results','dN_ptdptdphidy.dat')
for eta_s in eta_s_list:
	for tdec in tdec_list:
		# copy dNdyptdptdphi file from database to current iS if it exists
		dN_file_sourceFolder = path.join(backupDir_currentNode, 
			"etas_%g"%eta_s, "tdec_%g"%tdec)
		dN_file_sourceFile   = path.join(dN_file_sourceFolder,
			 "dN_ptdptdphidy.dat")
		if(path.isfile(dN_file_sourceFile)):
			shutil.rmtree(dN_file_iSfolder)# clear iS results folder
			mkdir(dN_file_iSfolder)
			for files in glob.glob(path.join(dN_file_sourceFolder, '*.*')):
				shutil.copy(files, dN_file_iSfolder)
			print "**** tdec = %g is calculating"%tdec
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
			print "**** tdec = %g finishs!"%tdec
		else:
			print "**** tdec = %g is skipped"%tdec
			continue
	print "runResDecay: eta/s=%g"%eta_s+" processed!"
print "All runs processed!"

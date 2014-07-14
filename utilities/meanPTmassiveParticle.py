#!/usr/bin/env python 

# file: meanPTmassiveParticle.py
# purpose: calculate the mean pT for a massive particle.
# how to use: keep it at the original folder of fs_package and run.

import numpy as np
from os import path, stat, getcwd
import sys

# run parameters, user define
particle_MCNum = 211 # Monte-Carlo number of hadron
matchingTime_list = np.linspace(0.4,1,4)
node_list = range(1,11)
events_per_node_list = range(1,41)

# system parameters
rootDir = getcwd() # save results to this folder
database_location = path.join(path.abspath('../..'),'dataBase') 
iStable_location = path.join(path.abspath('..'),'iS','tables')

def calculateMeanPT(particle_idx, pt_weight_array, targetFolder):
	"""
	calculate the mean PT for hadron with MC number particle_idx
	"""
	targetFile = path.join(targetFolder, "thermal_%d_vndata.dat"%particle_idx) # format:pt, m_T-m_0, dN/(dyptdpt2\pi)
	try:
		particle_result = np.loadtxt(targetFile)
	except:
		print "meanPTmassive: calculateMeanPT error!"
		print "Not such file as: "+targetFile
		sys.exit(-1)
	pt = particle_result[:,0]
	dNdyptdpt_2pi = particle_result[:,2]
	totalpt = 2.0*np.pi*sum(pt_weight_array*pt*pt*dNdyptdpt_2pi) #total pT=2pi*\int dpt*pt^2*dN/dyptdpt2pi
	totalnum = 2.0*np.pi*sum(pt_weight_array*pt*dNdyptdpt_2pi) # total num = 2pi*\int dpt*pt*dN/dyptdpt2pi
	result = totalpt/totalnum
	return(result)


def meanPTmassiveParticleShell():
	"""
	run over all nodes and specific matching times to find the meanPT of specific hadron
	"""
	# read in pt weight table for integration
	print "meanPTmassiveParticleShell: start!"

	table_file = path.join(iStable_location, 'pT_gauss_table.dat')
	pT_table = np.loadtxt(table_file) #format: pT, weight
	pT_weight = pT_table[:,1] 
	# construct a matrix for data 
	events_total = len(node_list)*len(events_per_node_list)
	mtimes_total = len(matchingTime_list)
	meanPT_table = np.zeros((events_total, mtimes_total)) 
	# start to loop over all events
	line_idx = 0
	for inode in node_list:
		for ievent in events_per_node_list:
			for itime_idx in range(mtimes_total):
				results_folder = path.join(database_location, 
					'node%d'%inode, 'event_%d'%ievent, '%g'%matchingTime_list[itime_idx])
				imeanPT = calculateMeanPT(particle_MCNum, pT_weight, results_folder)
				meanPT_table[line_idx, itime_idx] = imeanPT
			line_idx = line_idx+1 #next line for next event
	# save results table to a file for matlab plot
	meanPT_fileName = path.join(rootDir, 'meanPT_%d.dat'%particle_MCNum)
	commentLine = ''.join(["%10.2f"%i for i in matchingTime_list])
	np.savetxt(meanPT_fileName, meanPT_table, 
		fmt='%10.6e', comments='%',header=commentLine,
		delimiter='\t',) 
	print "meanPTmassiveParticleShell: finishes!"
	print "data save to file "+ meanPT_fileName


if __name__ == "__main__":
    meanPTmassiveParticleShell()



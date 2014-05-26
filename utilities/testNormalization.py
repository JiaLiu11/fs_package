#!/usr/bin/env python 
# Author: Jia Liu
# Purpose: analyze total energy at initial state and final state.
#
import dEcounters
import numpy as np
from os import path,stat
event_num = 99
matchingTimeList=np.linspace(8,10,3)
rootDir = path.abspath('..')  #define the root folder
table_location = path.join(rootDir, 'tables')
sfactor_list = np.loadtxt(path.join(table_location,'sfactor_log.dat'))

print 'tau', '	dEdy_thermal','	dEdy_hydro_init', '	dEdy_hydro_final',  '	dEdy_hydro_cf'
for tau in matchingTimeList:
	# calculate fs part total energy
	edFile = path.join(rootDir, 'fs','data', 'result', 'event_%d'%event_num, '%g'%tau,  'ed_profile_kln.dat')
	dEdyd2rdphipFolder = path.join(rootDir, 'fs','data', 'result', 'event_%d'%event_num, '%g'%tau)
	dEdyd2rdphipFile = path.join(dEdyd2rdphipFolder,'dEd2rdphip_kln.dat')
	sfactor = (sfactor_list[sfactor_list[:,0]==tau])[0,1]

	iSDataFolder = path.join(rootDir, 'localdataBase', 'event_%d'%event_num,'%g'%tau)
	targetFolder = iSDataFolder #back up data
	fileName='dEdydphipThermal.dat'
	dEdy_thermal = dEcounters.dEdydphipSigmaThermal(dEdyd2rdphipFile, edFile, sfactor, \
		targetFolder, 'dEdydphipThermal.dat')

	# calculate total initial energy
	dEdyphipFile = path.join(dEdyd2rdphipFolder,'dEdydphip_kln.dat')
	dEdy_init = dEcounters.dEdyTotalInital(dEdyphipFile, sfactor)

	# calculate final hydro energy
	iSFile = path.join(iSDataFolder,'decdat2.dat')
	if(stat(iSFile).st_size==0):
		dEdy_FO = 0
		dEdy_cf = 0
	else:
		dEdy_cf = dEcounters.dEdydphipSigmaFO('iS', iSDataFolder, targetFolder, 'dEdydphipFO.dat')
		dEdy_FO = dEcounters.dEdyTotalDecdat2(iSDataFolder)

	print tau, '\t', dEdy_thermal, '\t', dEdy_init-dEdy_thermal, '\t', dEdy_FO, '\t', dEdy_cf

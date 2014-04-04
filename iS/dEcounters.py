#!/usr/bin/env python 

#   file: runcode.py
#
#   Author:   Jia Liu    <liu.2053@osu.edu>
#   History:
#   Mar 03, 2014	 Include function getT0muFromiS() to calculate T0mu.
#   Feb 20, 2014	 Include rotation to T^0mu to convert 0 to local index
#   Feb 19, 2014     Disabled for now: Test mode: no initial flow: find 'jia test' for details.
#   Feb 18, 2014     Add dEdyViscousHydro() to calculate viscous hydro initial energy.
#   Feb 12, 2014     Add an argument iSFolder to dEdydphipSigmaFO(): let it find tables.
#   Feb 10, 2014     Add dEdyIdealHydro(): calculate initial ideal hydro part energy.
#   Feb 04, 2014     Add dEdyTotalInital(): calculate total initial energy.
#   Jan 31, 2014     dEdydphipSigmaFO() now use dE_dyptdptdphip table.
#   Jan 30, 2014     Add function dEdyTotalDecdat2() to calculate total energy
#				  flowing out of the freeze-out surface.
#   Jan 28, 2014     Fixed normalization bug of dEdydphipSigmaThermal.
#   Jan 24, 2014     Read the 'chosen_particles_backup.dat' to maintain
#    				 the correct matrix dimension. Commented dN calculation.

import sys
from os import path, stat
import numpy as np

# some constants
Edec = 0.18  #Unit: GeV/fm^3
dxdy = 0.01
includeBoost = True  # boost one index of T^munu to local rest frame

def  dEdydphipSigmaThermal(dEdyd2rdphipFile, edFile, sfactor, \
		targetFolder, fileName):
	"""
	kick out all the elements outside the freeze-out surface and do
	integration on the transverse plane. Return total energy of the
	thermalization surface and dEdydphip_{\Sigma_{thermal}}.
	"""
	dEdyd2rdphip_data = np.loadtxt(dEdyd2rdphipFile)
	ed_data = np.loadtxt(edFile)
	phipTbl = np.loadtxt('phip_gauss_table.dat') #100 points Gaussian points table for phip integration.
	phipGaussWeight = phipTbl[:, 1]

	# check dimensions
	ed_dataNum = ed_data.size
	if(ed_dataNum!=dEdyd2rdphip_data.shape[0]):
		print 'dEdydphipSigmaThermal error: Dimensions do not match!\n'
		print 'ed matrix elements: '+str(ed_dataNum), \
			', dEdyd2rdphip lines:' + str(dEdyd2rdphip_data.shape[0])
		sys.exit(-1)
	elif(phipGaussWeight.size!=dEdyd2rdphip_data.shape[1]):
		print 'dEdydphipSigmaThermal error: Dimensions do not match!\n'
		print 'phi_p Gaussian table length: '+str(phipGaussWeight.size),\
			', dEdyd2rdphip columns: ' + str(dEdyd2rdphip_data.shape[1])
		sys.exit(-1)

	# extract the row which are corresponding to cells outside the freeze-out 
	#surface
	ed_data = sfactor*ed_data
	ed_criteria = ed_data < Edec
	ed_criteria_rowNum = np.reshape(ed_criteria, (ed_dataNum)) #reshape is done row-wise
															# non-zero rows ed<Edec
	# integrate over the left cells
	dEdyd2rdphip_data = dEdyd2rdphip_data*sfactor
	dEdyd2rdphip_outside = dEdyd2rdphip_data[ed_criteria_rowNum, :]
	dEdydphip = dEdyd2rdphip_outside.sum(axis=0)*dxdy  # summing along each columnn
	dEdy = (dEdydphip*phipGaussWeight).sum()
	#print 'dEdyd2rdphip_outside', dEdyd2rdphip_outside.shape
	#print 'dEdydphip', dEdydphip.shape
	# save dEdydphip table
	savefileName = path.join(targetFolder, fileName)
	np.savetxt(savefileName, dEdydphip, fmt='%19.10e')
	return dEdy



def dEdydphipSigmaFO(iSFolder, iSDataFolder, targetFolder, fileName):
	"""
	calculate the total energy of the system from freeze-out surface (iS results)
	by summing over all thermal particles. Return the total transverse energy of
	freeze-out surface.
	"""
	# safty check, make sure iS has runned
	decdat2_file = path.join(iSDataFolder, 'decdat2.dat')
	if(stat(decdat2_file).st_size==0):  # empty decdat2 file, hydro does not run
		return 0	
	# get table folders and files
	EOSTblFolder = path.join(iSFolder, 'EOS')
	GaussTblFolder = path.join(iSFolder, 'tables')

	chosenParticleTbl = np.loadtxt(path.join(EOSTblFolder,'chosen_particles_backup.dat'))#debug
	ptTbl = np.loadtxt(path.join(GaussTblFolder, 'pT_gauss_table.dat'))
	phipTbl = np.loadtxt(path.join(GaussTblFolder, 'phi_gauss_table.dat'))
	dEdyptdptdphipTbl = np.loadtxt(path.join(iSDataFolder, 'dE_ptdptdphidy.dat'))
	particleMassTbl = np.loadtxt(path.join(EOSTblFolder,'particles_mass.dat'))
	particleMassTbl = particleMassTbl[:,2]  #first column is particle index

	#check the dimension of data tables
	particleTotalNum = chosenParticleTbl.size
	phipTblLength = phipTbl.shape[0] #debug
	ptTblLength = ptTbl.shape[0]
	# delete gamma from dNptdptdphidy table
	dEdyptdptdphipTbl=dEdyptdptdphipTbl[phipTblLength::,:]

	# dimension check
	if(dEdyptdptdphipTbl.shape[0]!=particleTotalNum*phipTblLength):
		print 'Wrong phi_p table!\n '
		print 'Length of phip Gaussian table:'+str(phipTblLength), \
			 ', dEdyptdptdphip table length: '+str(dEdyptdptdphipTbl.shape[0]/particleTotalNum) 
		sys.exit(-1)
	elif(dEdyptdptdphipTbl.shape[1]!=ptTblLength):
		print 'Wrong pT table: \n '
		print 'length of pT Gaussian table: '+str(ptTblLength),\
			', dEdyptdptdphip table length: '+str(dEdyptdptdphipTbl.shape[1])
		sys.exit(-1)

	dEdy = 0
	dEdydphipTbl = np.zeros((phipTblLength), float) #for all particles
	#dNdy = 0
	#dNdydphipTbl = np.zeros((phipTblLength), float) #for all particles
	pT_mat = np.tile(ptTbl[:,0].transpose(), (phipTblLength, 1))
	pTGaussWeight_mat = np.tile(ptTbl[:,1].transpose(), (phipTblLength, 1)) # tile(M, (m,n)): repeat matrix for m,n times
	phipGaussWeight_array = phipTbl[:,1]

	for i in range(0, particleTotalNum):
		tempTbl = dEdyptdptdphipTbl[i*phipTblLength:phipTblLength*(i+1), :] #per hadron
		# do integration on pT: dNTbl.*pT.*mT.*weight(pT)
		particleMass = particleMassTbl[i]
		dEdydphipTemp = (pT_mat*tempTbl*pTGaussWeight_mat).sum(axis=1)  #summing along one row
		# dNdydphipTemp = (pT_mat*tempTbl*pTGaussWeight_mat).sum(axis=1)
		#print 'Size of dEdydphip table:'+str(dEdydphipTemp.shape)
		#print 'Size of dEdydphip table:'+str(dEdydphipTbl.shape)
		dEdydphipTbl = dEdydphipTbl+dEdydphipTemp
		# dNdydphipTbl = dNdydphipTbl+dNdydphipTemp

		# integrate over phip
		dEdy = dEdy+(dEdydphipTemp*phipGaussWeight_array).sum()
		# dNdy = dNdy+(dNdydphipTemp*phipGaussWeight_array).sum()
		
	# save to file
	savefileName=path.join(targetFolder, fileName)
	np.savetxt(savefileName, dEdydphipTbl, fmt='%19.10e')

	# savefileName=path.join(targetFolder, 'dNdydphip.dat')
	# np.savetxt(savefileName, dNdydphipTbl, fmt='%19.10e')

	return dEdy


def dEdyTotalDecdat2(resultsFolder, fileName = 'decdat2.dat'):
	"""
	Calculate total energy from hydro freeze-out surface, i.e. file 'decdat2.dat'
	dEdy = \int T^{0\mu} d^3\sigma_\mu
	"""
	decdat2_file = path.join(resultsFolder, fileName)
	if(stat(decdat2_file).st_size==0):  # empty decdat2 file, hydro does not run
		return 0
	decData_fo = np.loadtxt(decdat2_file)

	# recombine T^0\mu
	vx_fo = decData_fo[:,4]
	vy_fo = decData_fo[:,5]
	gamma_fo = 1.0/np.sqrt(1.0 - vx_fo**2- vy_fo**2 )
	ux_fo = decData_fo[:,4]*gamma_fo
	uy_fo = decData_fo[:,5]*gamma_fo
	ed_fo = decData_fo[:,6]
	pl_fo = decData_fo[:,11]
	ppi_fo = decData_fo[:,19]
	pi00_fo = decData_fo[:,13]
	pi01_fo = decData_fo[:,14]
	pi02_fo = decData_fo[:,15]
	pi11_fo = decData_fo[:, 16]
	pi12_fo = decData_fo[:, 17]
	pi22_fo = decData_fo[:, 18]

	T00_fo = (ed_fo+pl_fo+ppi_fo)*gamma_fo*gamma_fo - pl_fo - ppi_fo \
		+ pi00_fo
	T01_fo = (ed_fo+pl_fo+ppi_fo)*gamma_fo*ux_fo + pi01_fo		
	T02_fo = (ed_fo+pl_fo+ppi_fo)*gamma_fo*uy_fo + pi02_fo 
	T0md3sigm = decData_fo[:,0]*(T00_fo*decData_fo[:,1] + T01_fo*decData_fo[:,2] \
		+ T02_fo*decData_fo[:,3])

	# include rotation? ! jia test
	global includeBoost
	if includeBoost == True:
		T11_fo = (ed_fo+pl_fo+ppi_fo)*ux_fo*ux_fo + pl_fo + ppi_fo \
			+ pi11_fo
		T12_fo = (ed_fo+pl_fo+ppi_fo)*ux_fo*uy_fo + pi12_fo
		T22_fo = (ed_fo+pl_fo+ppi_fo)*uy_fo*uy_fo + pl_fo + ppi_fo \
			+ pi22_fo
		T00_local = gamma_fo*T00_fo-ux_fo*T01_fo-uy_fo*T02_fo
		T01_local = gamma_fo*T01_fo-ux_fo*T11_fo-uy_fo*T12_fo
		T02_local = gamma_fo*T02_fo-ux_fo*T12_fo-uy_fo*T22_fo
		T0md3sigm = decData_fo[:,0]*(T00_local*decData_fo[:,1] + T01_local*decData_fo[:,2] \
			+ T02_local*decData_fo[:,3])

	dEdy = T0md3sigm.sum()

	return dEdy


def dEdytotalConstSurface(resultsFolder, fileName = 'constTauSurface.dat',dxdy=0.01):
	"""
	Calculate total energy from hydro constant time surface, i.e. file 'decdat2.dat'
	dEdeta_s = \int T^{00} d^3\sigma_0
	"""
	decdat2_file = path.join(resultsFolder, fileName)
	if(stat(decdat2_file).st_size==0):  # empty decdat2 file, hydro does not run
		return 0
	decData_fo = np.loadtxt(decdat2_file)

	# recombine T^0\mu
	vx_fo = decData_fo[:,4]
	vy_fo = decData_fo[:,5]
	gamma_fo = 1.0/np.sqrt(1.0 - vx_fo**2- vy_fo**2 )
	ux_fo = decData_fo[:,4]*gamma_fo
	uy_fo = decData_fo[:,5]*gamma_fo
	ed_fo = decData_fo[:,6]
	pl_fo = decData_fo[:,11]
	ppi_fo = decData_fo[:,19]
	pi00_fo = decData_fo[:,13]
	pi01_fo = decData_fo[:,14]
	pi02_fo = decData_fo[:,15]
	pi11_fo = decData_fo[:, 16]
	pi12_fo = decData_fo[:, 17]
	pi22_fo = decData_fo[:, 18]

	T00_fo = (ed_fo+pl_fo+ppi_fo)*gamma_fo*gamma_fo - pl_fo - ppi_fo \
		+ pi00_fo
	T00d3sigm = decData_fo[:,0]*T00_fo*dxdy
 
	# include boost?
	global includeBoost
	if includeBoost == True:
		T01_fo = (ed_fo+pl_fo+ppi_fo)*gamma_fo*ux_fo + pi01_fo		
		T02_fo = (ed_fo+pl_fo+ppi_fo)*gamma_fo*uy_fo + pi02_fo	
		T00_local = gamma_fo*T00_fo-ux_fo*T01_fo-uy_fo*T02_fo
		T00d3sigm = decData_fo[:,0]*T00_local*dxdy

	dEdy = T00d3sigm.sum()

	return dEdy




def dEdyTotalInital(dEdyphipFile, sfactor):
	"""
	calculate the total inital energy, by integrate the angular distribution dEdydphip
	"""
	try:
		dEdydphip_data = np.loadtxt(dEdyphipFile)
		phipTbl = np.loadtxt('phip_gauss_table.dat') #100 points Gaussian points table for phip integration.
	except:
		print 'No such files!'
		sys.exit(-1)
	phipGaussWeight = phipTbl[:, 1]
	dEdydphip_data*=sfactor
	# dimension check
	if(dEdydphip_data.shape[0]!=phipGaussWeight.shape[0]):
		print 'Wrong phi_p table!\n '
		print 'Length of phip Gaussian table:'+str(phipGaussWeight.shape[0]), \
			 ', dEdydphip table length: '+str(dEdydphip_data.shape[0]) 
		sys.exit(-1)

	dEdy = (dEdydphip_data*phipGaussWeight).sum()
	return dEdy


def dEdyIdealHydro(InitialFolder, sfactor, tau, dxdy):
	"""
	calculate the intial total energy inside hydro surface 
	if viscous terms are dropped directly.
	"""
	try:
		ed_data  = np.loadtxt(path.join(InitialFolder, 'ed_profile_kln.dat'))
		pressure_data = np.loadtxt(path.join(InitialFolder, 'Pressure_kln.dat'))
		ux_data = np.loadtxt(path.join(InitialFolder, 'ux_profile_kln.dat'))
		uy_data = np.loadtxt(path.join(InitialFolder, 'uy_profile_kln.dat'))
	except:
		print 'dEdyIdealHydro: reading files error!'
		sys.exit(-1)
	ed_scaled = ed_data*sfactor
	ed_criteria = ed_scaled >= Edec
	gamma_init = np.sqrt(1 + ux_data**2 + uy_data**2)
#	gamma_init = np.ones(gamma_init.shape) # jia test
	T00_init = (ed_data + pressure_data)*gamma_init*gamma_init \
		 - pressure_data
	T00_init_outside = T00_init[ed_criteria]  #only count out-of-freezeout elements
	dEdy_init = (T00_init_outside.sum()).sum()*dxdy*sfactor*tau
	return dEdy_init

def dEdyViscousHydro(InitialFolder, sfactor, tau, dxdy):
	"""
	calculate the intial total energy inside hydro surface 
	including viscous terms.
	"""
	try:
		ed_data  = np.loadtxt(path.join(InitialFolder, 'ed_profile_kln.dat'))
		pressure_data = np.loadtxt(path.join(InitialFolder, 'Pressure_kln.dat'))
		ux_data = np.loadtxt(path.join(InitialFolder, 'ux_profile_kln.dat'))
		uy_data = np.loadtxt(path.join(InitialFolder, 'uy_profile_kln.dat'))
		ppi_data = np.loadtxt(path.join(InitialFolder, 'BulkPi_kln.dat'))
		pi00_data = np.loadtxt(path.join(InitialFolder, 'Pi00_kln.dat'))
	except:
		print 'dEdyViscousHydro: reading files error!'
		sys.exit(-1)
	ed_scaled = ed_data*sfactor
	ed_criteria = ed_scaled >= Edec
	gamma_init = np.sqrt(1 + ux_data**2 + uy_data**2)
	T00_init = (ed_data + pressure_data + ppi_data)*gamma_init*gamma_init \
		 - pressure_data - ppi_data + pi00_data
	T00_init_outside = T00_init[ed_criteria]  #only count out-of-freezeout elements
	dEdy_init = (T00_init_outside.sum()).sum()*dxdy*sfactor*tau
	return dEdy_init

def getT0muFromiS(iSEOSFolder, iSDataFolder):
	"""
	Find T0mu of Cooper-Frye by summing over all thermal particles.
	"""
	try:
		chosenParticleTbl = np.loadtxt(path.join(iSEOSFolder, 'chosen_particles.dat'))
	except:
		print 'getT0muFromiS: cannot load chosen_particle table!'
	T0mu = np.zeros((4,1))
	i=0
	for particle_idx in chosenParticleTbl:
		thermalParticleETFile = 'thermal_%d_ET_integrated_vndata.dat' %particle_idx
		try:
			thermalParticleTbl = np.loadtxt(path.join(iSDataFolder, thermalParticleETFile))
		except:
			print 'getT0muFromiS: cannot load ET file for: '+str(particle_idx)
		for j in range(0,4):
			T0mu[j] += thermalParticleTbl[0, j+2]
		i+=1
	return T0mu[:]



#!/usr/bin/env python


# calculate total entropy 

import numpy as np
from scipy import interpolate
import os,sys

def loadTables(EOS_table_location, decdat2_location):
	"""
	Initialize the script by loading EOS table.
	Return a numpy array
	"""
	EOSTableFile = os.path.join(EOS_table_location, 'EOS_PST.dat')
	decdat2TableFile = os.path.join(decdat2_location, 'decdat2.dat')
	try:
		EOS_tbl = np.loadtxt(EOSTableFile)
	except:
		print 'loadTables: cannot load table:'+EOSTableFile
		sys.exit(-1)
	try:
		decdat2_tbl = np.loadtxt(decdat2TableFile)
	except:
		print 'loadTables: cannot load table:'+decdat2TableFile
		sys.exit(-1)		
	return EOS_tbl,decdat2_tbl


def convertEdToSd(ed_input, EOS_tbl):
	"""
	Convert energy density to entropy density.
	Input: ndarray, table value; 
	output: ndarray, interpolated.
	"""
	# extract the table which is used as the reference for interpolation
	ed_tbl = EOS_tbl[:, 0]
	sd_tbl = EOS_tbl[:, 2]
	# interpolate
	sd_output = np.interp(ed_input, ed_tbl, sd_tbl)
	return sd_output


def calculateTotalSFO(decdat2_tbl, EOS_tbl):
	"""
	calculate the total entropy flowing out of the freeze-out surface.
	"""
	# find entropy density
	ed_fo = decdat2_tbl[:, 6]
	sd_fo = convertEdToSd(ed_fo, EOS_tbl)
	# calculate entropy flow
	tauf = decdat2_tbl[:, 0]
	vx = decdat2_tbl[:, 4]
	vy = decdat2_tbl[:, 5]
	gamma = 1/np.sqrt(1-vx**2-vy**2)
	ux = gamma*vx
	uy = gamma*vy
	DA0 = decdat2_tbl[:, 1]
	DA1 = decdat2_tbl[:, 2]
	DA2 = decdat2_tbl[:, 3]
	umu_sigmamu = tauf*(DA0*gamma+DA1*ux+DA2*uy)
	stotal = (umu_sigmamu*sd_fo).sum()
	return stotal

def calculateTotoalSInit(hydro_init_folder, tau0, dxdy):
	"""
	calcualte initial total entropy.
	"""
	sd_init_file = os.path.join(hydro_init_folder, 'init-entropy.dat')
	sd_init = np.loadtxt(sd_init_file)
	sd_total = (tau0*dxdy*sd_init).sum()
	return sd_total


def totalEntrpyShell():
	EOS_tbl_loc = os.path.join('VISHNew','EOS','EOS_tables')
	decdat2_tbl_loc = os.path.join('iS','results')
	EOS_tbl, decdat2_tbl = loadTables(EOS_tbl_loc, decdat2_tbl_loc)

	# calculate total entropy on the freeze-out surface
	totalSd_FO = calculateTotalSFO(decdat2_tbl, EOS_tbl)
	print 'total entropy on freeze-out surface: %10.8e'%totalSd_FO
	totalSd_init = calculateTotoalSInit('VISHNew/Initial', 1.0, 0.01)
	print 'total entropy at initial:%10.8e'%totalSd_init


if __name__ == "__main__":
    totalEntrpyShell()

#!/usr/bin/env python 

from os import path, listdir, getcwd
import shutil

rootDir = getcwd()

def deleteAllFolders(folderName):
	Dir1 = [ x for x in listdir(rootDir) if path.isdir(x)]
	for i in Dir1:
		Dir2 = [x for x in listdir(path.join(rootDir, i)) if path.isdir(path.join(rootDir, i, x))]
		for j in Dir2:
			Dir3 = [x for x in listdir(path.join(rootDir, i, j)) if path.isdir(path.join(rootDir, i, j, x))]
			for k in Dir3:
				folderLocation = path.join(rootDir, i,j,k,\
					folderName)
				#print folderLocation
				if(path.isdir(folderLocation)):
					shutil.rmtree(folderLocation)
		print i, ' deleted!'


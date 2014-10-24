#!/usr/bin/env python

#   file: job_monitor.py
#
#   Author:   Jia Liu    <liu.2053@osu.edu>
#
#   This python script output the current progress of each job.

from os import path, getcwd, listdir

# some parameters
nodesList = range(1,11)

rootDir = getcwd();
dataBaseDir = path.join(rootDir, 'dataBase')

def getJobProcess(node_num):
    """
    return the latest event of one node.
    """	
    currentFolder = path.join(dataBaseDir, "node"+str(node_num))
    fileNames = listdir(currentFolder)
    try:
        fileNames.remove('Epx_initial.dat')  # eliminates additional names
    except:
        pass
    fileNames.remove('MCevents')
    fileNames = [items.replace('event_', '') for items in fileNames]
    fileNames = map(int, fileNames)
    return max(fileNames)


def getMatchingTime(node_num, event_num):
    """
    return the latest matching time of one node.
    """	
    currentFolder = path.join(dataBaseDir, "node"+ str(node_num), \
        'event_'+str(event_num));
    fileNames = listdir(currentFolder) 
    filesNames = map(float, fileNames)
    return max(fileNames)


def nodes_monitor_shell():
	"""
	Get the lastest finished event and matching time
	"""
	print 'Running report: '
	for i in nodesList:
		currentEvent = getJobProcess(i)
		currentMatchingTime = getMatchingTime(i, currentEvent)
		print 'Node ' + str(i) + ', event:' + str(currentEvent) \
			+ ' matching time: ' + str(currentMatchingTime)
        print 'Report finishes!\n'


if __name__ == "__main__":
    nodes_monitor_shell()

#! /usr/bin/env bash
# Put this file in the parent folder of main/

here=$(pwd)

if [ $# -lt 1 ];
then
    echo "Usage: runJobs.sh number_of_jobs"
    exit
fi	

for i in $(eval echo "{1..$1}"); do

	cd $here
	cd $here/nodes/node$i/utilities

	screen -dmS jia$i
	screen -S jia$i -p 0 -X stuff "python ./meanPT_totalN_calculator.py$(printf \\r)";
	
	echo "Nice work! Job $i submitted."
done;

echo "Done.";

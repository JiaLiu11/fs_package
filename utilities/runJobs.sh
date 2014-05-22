#! /usr/bin/env bash

here=$(pwd)

if [ $# -lt 1 ];
then
    echo "Usage: runJobs.sh number_of_jobs"
    exit
fi	

for i in $(eval echo "{1..$1}"); do

	cd $here
	rm -rf $here/nodes/node$i
	cp -r -a $here/main $here/nodes/node$i
	cd $here/nodes/node$i

	screen -dmS jia$i
	screen -S jia$i -p 0 -X stuff "python ./runcode.py$(printf \\r)";
	
	echo "Nice work! Job $i submitted."
done;

echo "Done.";

#! /usr/bin/env bash
(cd ../fs_package/
for ii in fs superMC VISHNew iS
    do
    (cd $ii; make; make clean)
done

echo "Compiling finished."
echo "Next generate jobs using generate-jobs-XXX.sh."
)

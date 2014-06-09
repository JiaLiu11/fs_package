#! /usr/bin/env bash
(cd ..
for ii in fs superMC VISHNew
    do
    (cd $ii; make; make clean)
done

echo "Compiling finished."
echo "Next generate jobs using generate-jobs-XXX.sh."
)

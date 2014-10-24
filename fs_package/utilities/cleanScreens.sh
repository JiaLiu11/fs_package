#! /usr/bin/env bash
# Author: Andy
# My style is to create a bunch of screens with my name, and then a descriptive word,
# for each job that I have to run.  This results in each receiving a number which will have
# the same count of digits as a prefix which represents each screen's unique id.  Using these ids,
# and providing the common tag words shared per job, you can quickly remove only the leftover
# screens from completed jobs.

if [ $# -lt 2 ]
then
    echo "Usage: cleanScreens.sh search_word number_of_digits_in_prefix"
    exit
fi

for i in $(screen -ls | grep $1 | cut -c 1-$(($2+1))); do screen -S $i -X quit; done


#!/bin/bash

file=$1
name=${file::-4}

echo "Residue,Occupancy" > ${name}.full.csv

# loop through numbers
for i in {1..521}
do
# check if number binds lipid (in file?)
        if grep -q "^$i," $1
        then
                grep "^$i," $1 >> ${name}.full.csv
        else
                echo ${i},0 >> ${name}.full.csv
        fi
# finish
done


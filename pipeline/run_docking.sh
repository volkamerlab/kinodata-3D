#!/bin/bash -l

# adjust to project directory
cd $HOME/kinodata-3D/

mkdir -p output/$1

start_time=$(date +%s.%N)

OE_LICENSE="./oe_license.txt" conda run --no-capture-output -n kinoml python pipeline/docking.py $1 data/templates/$2.pdb "$3" output/$1

end_time=$(date +%s.%N)

elapsed_time=$(( ${end_time%%.*} - ${start_time%%.*} ))


echo "Job terminated" >> output/$1/termination_reason.log
echo "runtime: $elapsed_time" >> output/$1/termination_reason.log


#!/bin/bash

for i in topol*; do
    awk '{ if (($1 == ";") && ($2 == "residue")) print $8}' $i > temp_${num}.txt
    num=$(($num + 1))
done

# cat the charges
cat temp_*.txt > charge_res.txt

# Current working directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Now we use a quick python script to sum the charges to give our required neutralisation charge
neu_charge=$(python $DIR/neutralise.py charge_res.txt)

# Clean up
rm temp_*.txt charge_res.txt

# Output charge
echo $neu_charge

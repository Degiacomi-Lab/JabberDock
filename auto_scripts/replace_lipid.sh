#!/bin/bash
# this should only run if the user has opted to used amber forcefield with slipid
loc=$(dirname $(echo "$0"))/memdat_files

# find and relace lipid17 resname with slipid
# packmol-memgen -l POPE:POPE -r 1:1 --salt --saltcon 0.0 --nocounter -p FILE.pdb --preoriented

input_file=${1}
output_file=${2}

# just in case file already exists
cp ${input_file} ${output_file}

# copy protein first
#${loc}/protein.awk ${input_file} >> ${output_file}
#echo "TER" >> ${output_file}

# change atomnames as appropiate
cat ${loc}/PA.dat | while read line; do
    convert_array=($line)
    ${loc}/lipid.awk ${output_file} $(printf '%4s' ${convert_array[0]}) $(printf '%4s' ${convert_array[1]}) "PA  " "POPE" > tmp.txt && mv tmp.txt ${output_file}
done

cat ${loc}/PE.dat | while read line; do
    convert_array=($line)
    ${loc}/lipid.awk ${output_file} $(printf '%4s' ${convert_array[0]}) $(printf '%4s' ${convert_array[1]}) "PE  " "POPE" > tmp.txt && mv tmp.txt ${output_file}
done
    
cat ${loc}/OL.dat | while read line; do
    convert_array=($line)
    ${loc}/lipid.awk ${output_file} $(printf '%4s' ${convert_array[0]}) $(printf '%4s' ${convert_array[1]}) "OL  " "POPE" > tmp.txt && mv tmp.txt ${output_file}
done

#fix water names
${loc}/water.awk ${output_file} > tmp.txt && mv tmp.txt ${output_file}

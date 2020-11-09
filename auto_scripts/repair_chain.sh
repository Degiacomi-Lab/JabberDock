#!/bin/bash

# get the fasta files needed to repair the pdbs listed in this folder
# then repair using autopatch 

loc=$(dirname $(echo "$0"))

# delimit by .
split_by () {
    string=$1
    separator=$2

    tmp=${string//"$separator"/$'\2'}
    IFS=$'\2' read -a arr <<< "$tmp"
    for substr in "${arr[@]}" ; do
        echo "$substr"
    done
    echo
}

struc=${1}
name=$(echo ${struc} | cut -c 1-4)
arr=($(split_by ${struc} '.'))
filename=$(echo ${arr[0]})

# need to only keep part of fasta that corresponds to our chain
if [[ ${#filename} =~ 4 ]]
then 
    chain=""
else
    chainarr=($(split_by ${filename} '_'))
    chain=$(echo ${chainarr[1]})
fi

if test -f ${filename}.fasta; then
    :
else
    wget https://www.rcsb.org/fasta/entry/${name}/download
    mv download ${filename}.fasta
fi

if [[ ${#chain} =~ 0 ]]
then 
    echo ""
else
    echo ">> Trimming fasta down to the specific chain(s)..."
    for (( i=0; i<${#chain}; i++ )) do
        letter="${chain:$i:1}"
	echo $letter
        awk -v C="Chain ${letter}" -v C2="Chains" '$0~C||$0~C2{print;getline;print}' ${filename}.fasta >> ${filename}_tmp.fasta
    done
    mv ${filename}_tmp.fasta ${filename}.fasta
fi

python ${loc}/autopatch.py ${filename} 15
python ${loc}/rmsd_align.py ${filename}_PATCHED.pdb ${filename}.pdb
mv ${filename}_PATCHED_align.pdb ${filename}.pdb

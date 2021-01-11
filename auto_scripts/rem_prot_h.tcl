mol new {fixed.pdb}

set sel [atomselect 0 "not (protein and hydrogen)"]
$sel writepdb fixed_trim.pdb

exit

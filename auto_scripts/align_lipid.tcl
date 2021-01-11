# TCL script to align trajectory along frame 1 (so no box jumping)
package require pbctools
mol new {npt_STID.gro}

set sel_gro [atomselect 0 "not resname SOL NA CL POPE"]
$sel_gro writepdb sel.pdb

mol delete 0

mol new {sel.pdb}
mol addfile {parse_traj.trr} type {trr} first 0 last -1 step 1 waitfor -1
animate delete beg 0 end 1 skip 0 1

# remove atom warping
pbc unwrap -sel "protein"
pbc wrap -compound fragment -center com -centersel "protein" -sel "protein" -all

# align along backbone
set start_sel [atomselect 1 "backbone" frame 0]
set current_sel [atomselect 1 "backbone"]
set all_sel [atomselect 1 "protein"]

set frames [molinfo 1 get numframes]
for { set f 0 } { $f < $frames } { incr f 1 } {

    $current_sel frame $f

    set trans_matrix [measure fit $current_sel $start_sel]
    
    $all_sel frame $f

    $all_sel move $trans_matrix
}

animate write pdb {./sim.pdb} beg 0 end -1 sel $all_sel

exit



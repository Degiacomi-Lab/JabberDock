#!/usr/bin/env python

""" A python script for writing a typical MD NVT / NPT input file for gromacs through user input specified in gmx_run.py, and the generate_ensemble input files used for dock.py """
import sys, os
import numpy as np

# Current path destination:
current_p = str(os.path.dirname(os.path.realpath(__file__)))

def nvtrun(ref_temp):
    line1 = "title   = NVT equib. simulation of protein via protocol program"
    line2 = "define   = -DPOSRES"
    line3 = "; Run Parameters"
    line4 = "integrator   = md"
    line5 = "nsteps   = 100000"
    line6 = "dt   =0.002"
    line7 = "; Output Control"
    line8 = "nstxout   = 1000"
    line9 = "nstvout   = 1000"
    line10 = "nstenergy   = 1000"
    line11 = "nstlog   = 1000"
    line12 = "; Bond Parameters"
    line13 = "continuation   = no"
    line14 = "constraint_algorithm   = lincs"
    line15 = "constraints   = all-bonds"
    line16 = "lincs_iter   = 1"
    line17 = "lincs_order   = 4"
    line18 = "; Neighbour Searching"
    line18_5 = "cutoff-scheme   = Verlet"
    line19 = "ns_type   = grid"
    line20 = "nstlist   = 5"
    line21 = "rlist   = 1.2 ; Cutoff range for neighbourlist (nm)"
    line22 = "rcoulomb   = 1.2 ; Cutoff range for electrostatics (nm)"
    line23 = "rvdw   = 1.2 ; Cutoff range for vdW (nm)"
    line24 = "; Electrostatics"
    line25 = "coulombtype   = PME"
    line26 = "pme_order   = 4"
    line27 = "fourierspacing   = 0.16"
    line28 = "; Temp coulpling is on"
    line29 = "tcoupl   = V-rescale"
    line30 = "tc-grps   = Protein Non-Protein" # This one needs checking for salt
    line31 = "tau_t   = 0.1 0.1" # Coupling constants between above groups
    line32 = "ref_t   = %f %f"%(ref_temp, ref_temp)
    line33 = "; Pressure coupling is off"
    line34 = "pcoupl   = no"
    line39 = "; PBC"
    line40 = "pbc   = xyz"
    line41 = "; Dispersion Correction"
    line42 = "DispCorr   = EnerPres"
    line43 = "; Velocity Generation"
    line44 = "gen_vel   = yes"
    line45 = "gen_temp   = %f"%(ref_temp)
    line46 = "gen_seed   = -1"

    lines = [line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12, line13, line14, line15, line16, line17, line18, line18_5, line19, line20, line21, line22, line23, line24, line25, line26, line27, line28, line29, line30, line31, line32, line33, line34, line39, line40, line41, line42, line43, line44, line45, line46]

    md_file = open('nvt.mdp', 'w')

    for line in lines:
        line += '\n'
        md_file.write(line)
    md_file.close()

    return 0

def nptrun(timestep, no_steps, dump_time, ref_temp, ref_press):

    timestep /= 1000.

    line1 = "title   = NPT simulation of protein via protocol program"
    line2 = "; No posres used as we're running from an equib sim."
    line3 = "; Run Parameters"
    line4 = "integrator   = md"
    line5 = "nsteps   = %i"%(no_steps)
    line6 = "dt   = %f"%(timestep)
    line7 = "; Output Control"
    line8 = "nstxout   = %i"%(dump_time)
    line9 = "nstvout   = %i"%(dump_time)
    line10 = "nstenergy   = %i"%(dump_time)
    line11 = "nstlog   = %i"%(dump_time)
    line12 = "; Bond Parameters"
    line13 = "continuation   = yes"
    line14 = "constraint_algorithm   = lincs"
    line15 = "constraints   = all-bonds"
    line16 = "lincs_iter   = 1"
    line17 = "lincs_order   = 4"
    line18_5 = "cutoff-scheme   = Verlet"
    line18 = "; Neighbour Searching"
    line19 = "ns_type   = grid"
    line20 = "nstlist   = 5"
    line21 = "rlist   = 1.2 ; Cutoff range for neighbourlist (nm)"
    line22 = "rcoulomb   = 1.2 ; Cutoff range for electrostatics (nm)"
    line23 = "rvdw   = 1.2 ; Cutoff range for vdW (nm)"
    line24 = "; Electrostatics"
    line25 = "coulombtype   = PME"
    line26 = "pme_order   = 4"
    line27 = "fourierspacing   = 0.16"
    line28 = "; Temp coulpling is on"
    line29 = "tcoupl   = Berendsen"
    line30 = "tc-grps   = Protein Non-Protein" # This one needs checking for salt
    line31 = "tau_t   = 0.1 0.1" # Coupling constants between above groups
    line32 = "ref_t   = %f %f"%(ref_temp, ref_temp)
    line33 = "; Pressure coupling is on"
    line34 = "pcoupl   = Berendsen" # Simple barostat, don't really need time dynamics
    line35 = "pcoupltype   = semiisotropic" # isotropic can be used for non-membrane proteins
    line36 = "tau_p   = 10.0"
    line37 = "ref_p   = %f %f"%(ref_press, ref_press)
    line38 = "compressibility = 4.5e-5 4.5e-5"
    line39 = "; PBC"
    line40 = "pbc   = xyz"
    line41 = "; Dispersion Correction"
    line42 = "DispCorr   = EnerPres"
    line43 = "; Velocity Generation"
    line44 = "gen_vel   = no"
    #line45 = "gen_temp   = %f"%(ref_temp)
    #line46 = "gen_seed   = -1"
    line47 = "; COM motion removal, left intentionally commented out as are used with membrane proteins"
    line48 = "; nstcomm = 1"
    line49 = "; comm-mode = linear"
    line50 = "; comm-grps = Protein"
    line51 = "; Scale COM of reference coordinates"
    line52 = "; refcoord_scaling = com"

    lines = [line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12, line13, line14, line15, line16, line17, line18, line18_5, line19, line20, line21, line22, line23, line24, line25, line26, line27, line28, line29, line30, line31, line32, line33, line34, line35, line36, line37, line38, line39, line40, line41, line42, line43, line44, line47, line48, line49, line50, line51, line52]

    md_file = open('npt.mdp', 'w')

    for line in lines:
        line += '\n'
        md_file.write(line)

    md_file.close()
    return 0

def nptrun_mem(timestep, no_steps, dump_time, ref_temp, ref_press, continuation="no", ions=True, lipid="POPE", constrain=False, equib=True):

    timestep /= 1000.

    line1 = "title   = NPT simulation of protein in membrane via protocol program"
    if constrain:
        line2 = "define  = -DPOSRES  ; position restrain the membrane"
    else:
        line2 = ";define  = -DPOSRES  ; position restrain the membrane"
    line3 = "; Run Parameters"
    line4 = "integrator   = md"
    line5 = "nsteps   = %i"%(no_steps)
    line6 = "dt   = %f"%(timestep)
    line7 = "; Output Control"
    line8 = "nstxout   = %i"%(dump_time)
    line9 = "nstvout   = %i"%(dump_time)
    line10 = "nstenergy   = %i"%(dump_time)
    line11 = "nstlog   = %i"%(dump_time)
    line12 = "; Bond Parameters"
    line13 = "continuation   = %s"%(continuation)
    line14 = "constraint_algorithm   = lincs"
    line15 = "constraints   = h-bonds"
    line16 = "lincs_iter   = 1"
    line17 = "lincs_order   = 4"
    line18_5 = "cutoff-scheme   = Verlet"
    line18 = "; Neighbour Searching"
    line19 = "ns_type   = grid"
    line20 = "nstlist   = 5"
    line21 = "rlist   = 1.2 ; Cutoff range for neighbourlist (nm)"
    line22 = "rcoulomb   = 1.2 ; Cutoff range for electrostatics (nm)"
    line23 = "rvdw   = 1.2 ; Cutoff range for vdW (nm)"
    line24 = "; Electrostatics"
    line25 = "coulombtype   = PME"
    line26 = "pme_order   = 4"
    line27 = "fourierspacing   = 0.16"
    line28 = "; Temp coulpling is on"
    line29 = "tcoupl   = V-rescale"
    if ions:
        line30 = "tc-grps   = protein %s Water_and_ions"%(lipid) # This one needs checking for salt
    else:
        line30 = "tc-grps   = protein %s Water"%(lipid) 
    line31 = "tau_t   = 0.1 0.1 0.1" # Coupling constants between above groups
    line32 = "ref_t   = %f %f %f"%(ref_temp, ref_temp, ref_temp)
    line33 = "; Pressure coupling is on"
    line34 = "pcoupl   = Berendsen" # Simple barostat, don't really need time dynamics
    line35 = "pcoupltype   = semiisotropic" # isotropic can be used for non-membrane proteins
    line36 = "tau_p   = 1.0"
    line37 = "ref_p   = %f %f"%(ref_press, ref_press)
    line38 = "compressibility = 4.5e-5 4.5e-5"
    line38_5 = "refcoord_scaling        = com"
    line39 = "; PBC"
    line40 = "pbc   = xyz"
    line41 = "; Dispersion Correction"
    line42 = "DispCorr   = EnerPres"
    line43 = "; Velocity Generation"
    if continuation == "no":
        line44 = "gen_vel = yes"
        line45 = "gen_temp  = %f   ; temperature for Maxwell distribution"%(ref_temp)
        line46 = "gen_seed  = -1    ; generate a random seed"
    elif continuation == "yes":
        line44 = "gen_vel = no"
        line45 = ""
        line46 = ""

    lines = [line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12, line13, line14, line15, line16, line17, line18, line18_5, line19, line20, line21, line22, line23, line24, line25, line26, line27, line28, line29, line30, line31, line32, line33, line34, line35, line36, line37, line38, line38_5, line39, line40, line41, line42, line43, line44, line46]

    if equib:
        if constrain:
            md_file = open('npt_constrain.mdp', 'w')
        else:
            md_file = open('npt_noconstrain.mdp', 'w')
    else:
        md_file = open('npt_STID.mdp', 'w')

    for line in lines:
        line += '\n'
        md_file.write(line)

    md_file.close()
    return 0

def powrun(x, y, z, receptor, ligand, iso = 0.43, dist = 1.6, log='pow_log.dat',fname='input_ensemble', no_samples = 300, angle=180., mpi=False, axis=1, file_name = 'model_solutions', restart=False, ensemble_module="generate_ensemble.py"):
    
    # axis refers to limiting the rotational space of protein (needs to be added in docs)
    # ensemble module is the generate ensemble pow script that someone can edit (needs to be added in docs)

    receptor_pdb = receptor + '.pdb'
    ligand_pdb = ligand + '.pdb'
    receptor_dx = receptor + '.dx'
    ligand_dx = ligand + '.dx'

    line1 = "optimizer PSO"
    line2 = "module %s/%s"%(current_p, ensemble_module)
    line3 = " "
    line4 = "# optimizer parameters"
    line5 = "steps 300"
    line6 = "particles 80"
    line7 = "repeat 3"
    line8 = "repulsion on"
    line9 = "kar_threshold 0.01"
    line10 = "repulsion_factor 0.04"
    line11 = " "
    line12 = "# definition of search space boundary conditions (crd, angle, axis of rotation)"
    line13 = "boundaryMin -%f -%f -%f -%f -%f -%f -%f"%(x, y, z, angle, axis, axis, axis)
    line14 = "boundaryMax %f %f %f %f %f %f %f"%(x, y, z, angle, axis, axis, axis)
    line15 = " "
    line16 = "# name of logfile"
    line17 = "output %s"%(log)
    line18_5 = " "
    line18 = "# Filtering criterion from data in logfile"
    line19 = "filter_threshold 0.0"
    line20 = " "
    line21 = "# CUSTOM KEYWORDS"
    line22 = "receptor %s"%(receptor_pdb)
    line23 = "ligand %s"%(ligand_pdb)
    line24 = "samples %i"%(no_samples)
    line26 = "file_name %s"%(file_name)
    line27 = "receptor_map %s" %(receptor_dx)
    line28 = "ligand_map %s"%(ligand_dx)
    line29 = "iso_cutoff %s"%(iso)
    line30 = "dist_cutoff %s"%(dist)
    if restart:
        line31 = "restart swarm.restart"
    else:
        line31 = ""

    lines = [line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12, line13, line14, line15, line16, line17, line18, line18_5, line19, line20, line21, line22, line23, line24, line26, line27, line28, line29, line30, line31]

    md_file = open('%s'%(fname), 'w')

    for line in lines:
        line += '\n'
        md_file.write(line)

    md_file.close()

    return 0

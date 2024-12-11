#!/bin/bash
#
# Which partition to use
#MSUB -q skylake
# Number of nodes
#MSUB -n 96
#MSUB -c 2
# Time limit in seconds
#MSUB -T 86000
# Job name
#MSUB -r slurm
# Mail alerts
#MSUB -@ j.h.kasper@utwente.nl:begin,end
#

#ccc_mprun LESsolver
ccc_mprun -f mpmd_irene.conf

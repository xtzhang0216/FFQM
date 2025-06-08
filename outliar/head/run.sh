#!/bin/bash
#SBATCH --output /pubhome/xtzhang/interaction/FF/outliar/outliar_head/MIMM_N1PA/939_%j.out
#SBATCH --error /pubhome/xtzhang/interaction/FF/outliar/outliar_head/MIMM_N1PA/939_%j.err
#SBATCH --job-name qm
#SBATCH --partition benz
#SBATCH --nodes 1
#SBATCH --ntasks 32
#SBATCH --mem 250G
#SBATCH --nodelist k154
#SBATCH --qos long_few
#SBATCH --time 3-00:00:00
source ~/.bashrc
cd /tmp/xtzhang
test -d MIMM_N1PA || mkdir MIMM_N1PA
cd MIMM_N1PA
cp /pubhome/xtzhang/interaction/FF/outliar/outliar_head/MIMM_N1PA/939.inp .
/pubhome/soft/orca/orca_6_0_0/orca 939.inp > 939.log
cp 939.log /pubhome/xtzhang/interaction/FF/outliar/outliar_head/MIMM_N1PA
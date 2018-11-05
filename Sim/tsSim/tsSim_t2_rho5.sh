#!/bin/bash
#$ -o /storage/essicd/data/HCP/Soroosh/ACAnal/tsSim/logs
#$ -e /storage/essicd/data/HCP/Soroosh/ACAnal/tsSim/logs
#$ -l h_rt=00:10:00
#$ -l h_vmem=3G
#$ -N tsSim_t2_rho5
#$ -r y
#$ -t 1-1000

cd /home/wmrnaq/ACAnal/Sim/tsSim

echo "t2_rho5"

.  /etc/profile
module add matlab
matlab -nodisplay -nosplash -nodesktop -nojvm -r "rho_cnt=5;t_cnt=2;tsSim"

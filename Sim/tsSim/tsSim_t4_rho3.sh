#!/bin/bash
#$ -o /storage/essicd/data/HCP/Soroosh/ACAnal/tsSim/logs
#$ -e /storage/essicd/data/HCP/Soroosh/ACAnal/tsSim/logs
#$ -l h_rt=00:10:00
#$ -l h_vmem=3G
#$ -N tsSim_t4_rho3
#$ -r y
#$ -t 774

cd /home/wmrnaq/ACAnal/Sim/tsSim

echo "t4_rho3"

.  /etc/profile
module add matlab
matlab -nodisplay -nosplash -nodesktop -nojvm -r "rho_cnt=3;t_cnt=4;tsSim"

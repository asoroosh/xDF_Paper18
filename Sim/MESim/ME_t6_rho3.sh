#!/bin/bash
#$ -o /storage/essicd/data/HCP/Soroosh/ACAnal/MESim/logs
#$ -e /storage/essicd/data/HCP/Soroosh/ACAnal/MESim/logs
#$ -q short.q
#$ -l h_rt=00:35:00
#$ -l h_vmem=5G
#$ -N MEq_t6_rho3
#$ -r y
#$ -t 1-2000

cd /home/wmrnaq/ACAnal/Sim/MESim

.  /etc/profile

module add matlab
matlab -nodisplay -nosplash -nodesktop -nojvm -r "rho_cnt=3;t_cnt=6;ME_Sim"

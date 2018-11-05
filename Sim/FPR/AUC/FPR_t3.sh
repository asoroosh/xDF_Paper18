#!/bin/bash
#$ -o /storage/essicd/data/HCP/Soroosh/ACAnal/MESim/FPR/SenSpc/logs
#$ -e /storage/essicd/data/HCP/Soroosh/ACAnal/MESim/FPR/SenSpc/logs
#$ -l h_rt=02:00:00
#$ -l h_vmem=6G
#$ -N SenSpc_AUC_3
#$ -r y
#$ -t 1-2000

cd /home/wmrnaq/ACAnal/Sim/FPR/AUC

.  /etc/profile

module add matlab
matlab -nodisplay -nosplash -nodesktop -nojvm -r "t_cnt=3;SenSpc"

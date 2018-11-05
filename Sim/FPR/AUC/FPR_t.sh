#!/bin/bash
#$ -o /home/wmrnaq/ACAnal/Sim/FPR/AUC/logs
#$ -e /home/wmrnaq/ACAnal/Sim/FPR/AUC/logs
#$ -l h_rt=01:00:00
#$ -l h_vmem=6G
#$ -N SenSpc_4
#$ -r y
#$ -t 1

cd /home/wmrnaq/ACAnal/Sim/FPR/AUC

.  /etc/profile

module add matlab
matlab -nodisplay -nosplash -nodesktop -nojvm -r "t_cnt=4;SenSpc"

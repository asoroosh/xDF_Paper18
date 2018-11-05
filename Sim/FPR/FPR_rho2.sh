#!/bin/bash
#$ -o /home/wmrnaq/ACAnal/Sim/FPR/logs
#$ -e /home/wmrnaq/ACAnal/Sim/FPR/logs
#$ -l h_rt=01:00:00
#$ -l h_vmem=6G
#$ -N SenSpc_2
#$ -r y
#$ -t 1-2000

cd /home/wmrnaq/ACAnal/Sim/FPR

.  /etc/profile

module add matlab
matlab -nodisplay -nosplash -nodesktop -nojvm -r "rho_cnt=2;FPR_VarRho"

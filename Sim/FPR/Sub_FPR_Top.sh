#!/bin/bash
#$ -o /home/wmrnaq/ACAnal/Sim/FPR/logs
#$ -e /home/wmrnaq/ACAnal/Sim/FPR/logs
#$ -q short.q
#$ -l h_rt=01:00:00
#$ -l h_vmem=6G
#$ -N SenSpc
#$ -r y
#$ -t 1-1000

cd /home/wmrnaq/ACAnal/Sim/FPR

.  /etc/profile
module add matlab

matlab -nodisplay -nojvm -nosplash -nodesktop -r "FPR_Top"

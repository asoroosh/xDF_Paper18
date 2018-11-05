#!/bin/bash
#$ -o /home/wmrnaq/ACAnal/Sim/FPR/AUC/
#$ -e /home/wmrnaq/ACAnal/Sim/FPR/AUC/
#$ -l h_rt=48:00:00
#$ -l h_vmem=12G
#$ -N Agg_AUC
#$ -r y

cd /home/wmrnaq/ACAnal/Sim/FPR/AUC

.  /etc/profile

module add matlab
matlab -nodisplay -nosplash -nodesktop -nojvm -r "Agg_AUC"

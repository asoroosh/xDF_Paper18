#!/bin/bash
#$ -o /home/wmrnaq/ACAnal/MassAC/logs
#$ -e /home/wmrnaq/ACAnal/MassAC/logs
#$ -l h_rt=02:00:00
#$ -l h_vmem=20G
#$ -N MassAC
#$ -r y
#$ -t 1-100

cd /home/wmrnaq/ACAnal/MassAC

.  /etc/profile
module add matlab

matlab -nodisplay -nosplash -nodesktop -r "HCP_100UR_MassAC"

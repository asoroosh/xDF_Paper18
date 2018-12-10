#!/bin/bash
#$ -o /home/wmrnaq/HCP/FC/Mats/DenRng/logs
#$ -e /home/wmrnaq/HCP/FC/Mats/DenRng/logs
#$ -l h_rt=00:40:00
#$ -l h_vmem=4G
#$ -N FCGM_Den8_Gordon_FPP_100UR
#$ -r y
#$ -t 1-100

cd ~/HCP/FC/Mats/DenRng

.  /etc/profile
module add matlab

matlab -nojvm -nodisplay -nosplash -nodesktop -r "PP='FPP';DenID=8;Atlas='Gordon';FCGMDenRng;"

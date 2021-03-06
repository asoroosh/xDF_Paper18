#!/bin/bash
#$ -o /home/wmrnaq/HCP/FC/Mats/CE/logs
#$ -e /home/wmrnaq/HCP/FC/Mats/CE/logs
#$ -l h_rt=00:40:00
#$ -l h_vmem=4G
#$ -N Gordon_MPP_GM_100UR
#$ -r y
#$ -t 1-100

cd ~/HCP/FC/Mats/CE

.  /etc/profile
module add matlab

matlab -nojvm -nodisplay -nosplash -nodesktop -r "PP='MPP';Atlas='Gordon';FCGMCE;"

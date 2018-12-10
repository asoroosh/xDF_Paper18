#!/bin/bash
#$ -o /home/wmrnaq/HCP/FC/Mats/Stat/logs
#$ -e /home/wmrnaq/HCP/FC/Mats/Stat/logs
#$ -l h_rt=00:40:00
#$ -l h_vmem=5G
#$ -N Yeo_FPP_GM_100UR
#$ -r y
#$ -t 1-100

cd ~/HCP/FC/Mats/Stat

.  /etc/profile
module add matlab

matlab -nojvm -nodisplay -nosplash -nodesktop -r "PP='FPP';Atlas='Yeo';FCGMStat;"

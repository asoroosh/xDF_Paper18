#!/bin/bash
#$ -o /home/wmrnaq/HCP/FC/Mats/Unthr/logs
#$ -e /home/wmrnaq/HCP/FC/Mats/Unthr/logs
#$ -l h_rt=00:60:00
#$ -l h_vmem=4G
#$ -N Yeo_FPP_GM_100UR
#$ -r y
#$ -t 1-100

cd ~/HCP/FC/Mats/Unthr

.  /etc/profile
module add matlab

matlab -nojvm -nodisplay -nosplash -nodesktop -r "PP='FPP';Atlas='Yeo';FCGMUnthr;"

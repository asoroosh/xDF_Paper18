#!/bin/bash
#$ -o /home/wmrnaq/DyConn/InterSubMats/logs
#$ -e /home/wmrnaq/DyConn/InterSubMats/logs
#$ -l h_rt=01:30:00
#$ -l h_vmem=5G
#$ -N YEO_DyConnInterSub_100UR
#$ -r y
#$ -t 1-5000

cd /home/wmrnaq/DyConn/InterSubMats

.  /etc/profile
module add matlab

matlab -nojvm -nodisplay -nosplash -nodesktop -r "Atlas='YEO';InterSubMat"

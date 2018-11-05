#!/bin/bash
#$ -o /home/wmrnaq/ACAnal/Sim/Topology/logs
#$ -e /home/wmrnaq/ACAnal/Sim/Topology/logs
#$ -q short.q
#$ -l h_rt=01:00:00
#$ -l h_vmem=12G
#$ -N Yeo_Null_Mats
#$ -r y
#$ -t 1-2000

cd /home/wmrnaq/ACAnal/Sim/Topology

.  /etc/profile
module add matlab

matlab -nodisplay -nosplash -nodesktop -r "Yeo_2_FakeMats"

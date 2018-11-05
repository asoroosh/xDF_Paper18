
#!/bin/bash
#$ -o /home/wmrnaq/ACAnal/Sim/MESim/Roy/logs
#$ -e /home/wmrnaq/ACAnal/Sim/MESim/Roy/logs
#$ -l h_rt=00:20:00
#$ -l h_vmem=5G
#$ -N RoyRoch_1_4
#$ -r y
#$ -t 1-2000

.  /etc/profile
module load R/3.3.2
cd /home/wmrnaq/ACAnal/Sim/MESim/Roy
/usr/local/packages/R-3.3.2/bin/R --vanilla --restore  t_cnt=1 r_cnt=4 < Sim_RoyVarEst.R

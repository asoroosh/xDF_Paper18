
for t_cnt in `seq 1 6`
do
  	for rho_cnt in `seq 1 5`
        do
FileName="Roy_t${t_cnt}_rho${rho_cnt}.sh"
cat > $FileName << EOF

#!/bin/bash
#$ -o /home/wmrnaq/ACAnal/Sim/MESim/Roy/logs
#$ -e /home/wmrnaq/ACAnal/Sim/MESim/Roy/logs
#$ -l h_rt=00:20:00
#$ -l h_vmem=5G
#$ -N RoyRoch_${t_cnt}_${rho_cnt}
#$ -r y
#$ -t 1-2000

.  /etc/profile
module load R/3.3.2
cd /home/wmrnaq/ACAnal/Sim/MESim/Roy
/usr/local/packages/R-3.3.2/bin/R --vanilla --restore  t_cnt=$t_cnt r_cnt=$rho_cnt < Sim_RoyVarEst.R

EOF
                qsub $FileName
        done
done






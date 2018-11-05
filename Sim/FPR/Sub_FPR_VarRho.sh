  	for rho_cnt in `seq 1 5`
        do
FileName="FPR_rho${rho_cnt}.sh"
cat > $FileName << EOF
#!/bin/bash
#$ -o /home/wmrnaq/ACAnal/Sim/FPR/logs
#$ -e /home/wmrnaq/ACAnal/Sim/FPR/logs
#$ -l h_rt=01:00:00
#$ -l h_vmem=6G
#$ -N SenSpc_${rho_cnt}
#$ -r y
#$ -t 1-2000

cd /home/wmrnaq/ACAnal/Sim/FPR

.  /etc/profile

module add matlab
matlab -nodisplay -nosplash -nodesktop -nojvm -r "rho_cnt=${rho_cnt};FPR_VarRho"
EOF
                qsub $FileName
        done


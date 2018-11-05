  	for t_cnt in `seq 1 2`
        do
FileName="FPR_t${t_cnt}.sh"
cat > $FileName << EOF
#!/bin/bash
#$ -o /home/wmrnaq/ACAnal/Sim/FPR/AUC/logs
#$ -e /home/wmrnaq/ACAnal/Sim/FPR/AUC/logs
#$ -l h_rt=03:00:00
#$ -l h_vmem=20G
#$ -N SenSpc_AUC_${t_cnt}
#$ -r y
#$ -t 1-2000

cd /home/wmrnaq/ACAnal/Sim/FPR/AUC

.  /etc/profile

module add matlab
matlab -nodisplay -nosplash -nodesktop -nojvm -r "t_cnt=${t_cnt};SenSpc_LongTS"
EOF
                qsub $FileName
        done


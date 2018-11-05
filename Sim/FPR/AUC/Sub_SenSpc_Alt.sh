  	for t_cnt in `seq 1 4`
        do
FileName="FPR_t${t_cnt}_Alt.sh"
cat > $FileName << EOF
#!/bin/bash
#$ -o /home/wmrnaq/ACAnal/Sim/FPR/AUC/logs
#$ -e /home/wmrnaq/ACAnal/Sim/FPR/AUC/logs
#$ -l h_rt=02:00:00
#$ -l h_vmem=6G
#$ -N SenSpc_AUC_${t_cnt}_Alt
#$ -r y
#$ -t 1-500

cd /home/wmrnaq/ACAnal/Sim/FPR/AUC

.  /etc/profile

module add matlab
matlab -nodisplay -nosplash -nodesktop -nojvm -r "t_cnt=${t_cnt};SenSpc_Alt"
EOF
                qsub $FileName
        done


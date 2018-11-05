  	for t_cnt in `seq 1 4`
        do
FileName="FPR_t${t_cnt}.sh"
cat > $FileName << EOF
#!/bin/bash
#$ -o /storage/essicd/data/HCP/Soroosh/ACAnal/MESim/FPR/SenSpc/logs
#$ -e /storage/essicd/data/HCP/Soroosh/ACAnal/MESim/FPR/SenSpc/logs
#$ -l h_rt=02:00:00
#$ -l h_vmem=6G
#$ -N SenSpc_AUC_${t_cnt}
#$ -r y
#$ -t 1-2000

cd /home/wmrnaq/ACAnal/Sim/FPR/AUC

.  /etc/profile

module add matlab
matlab -nodisplay -nosplash -nodesktop -nojvm -r "t_cnt=${t_cnt};SenSpc"
EOF
                qsub $FileName
        done


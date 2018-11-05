nRlz=5000

mat_filename="ME_Sim_TVOff"

for t_cnt in `seq 1 6`
do
	for rho_cnt in `seq 1 5`
	do
FileName="ME_t${t_cnt}_rho${rho_cnt}.sh"
cat > $FileName << EOF
#!/bin/bash
#$ -o /storage/essicd/data/HCP/Soroosh/ACAnal/MESim/logs
#$ -e /storage/essicd/data/HCP/Soroosh/ACAnal/MESim/logs
#$ -q short.q
#$ -l h_rt=00:35:00
#$ -l h_vmem=5G
#$ -N MEq_t${t_cnt}_rho${rho_cnt}
#$ -r y
#$ -t 1-${nRlz}

cd /home/wmrnaq/ACAnal/Sim/MESim

.  /etc/profile


module add matlab
matlab -nodisplay -nosplash -nodesktop -nojvm -r "rho_cnt=${rho_cnt};t_cnt=${t_cnt};${mat_filename}"
EOF
		qsub $FileName
        done
done

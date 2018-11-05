nRlz=1000

mat_filename="tsSim"

for t_cnt in `seq 1 4`
do	
	for rho_cnt in `seq 1 5`
	do
FileName="tsSim_t${t_cnt}_rho${rho_cnt}.sh"
cat > $FileName << EOF
#!/bin/bash
#$ -o /storage/essicd/data/HCP/Soroosh/ACAnal/tsSim/logs
#$ -e /storage/essicd/data/HCP/Soroosh/ACAnal/tsSim/logs
#$ -l h_rt=00:10:00
#$ -l h_vmem=3G
#$ -N tsSim_t${t_cnt}_rho${rho_cnt}
#$ -r y
#$ -t 1-${nRlz}

cd /home/wmrnaq/ACAnal/Sim/tsSim

echo "t${t_cnt}_rho${rho_cnt}"

.  /etc/profile
module add matlab
matlab -nodisplay -nosplash -nodesktop -nojvm -r "rho_cnt=${rho_cnt};t_cnt=${t_cnt};${mat_filename}"
EOF

		qsub $FileName
        done
done

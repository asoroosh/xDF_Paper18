nRlz=2000

mat_filename="FPR_ME_Sim"

for t_cnt in `seq 1 4`
do
	for rho_cnt  in  1
	do
FileName="ME_t${t_cnt}_rho${rho_cnt}.sh"
cat > $FileName << EOF
#!/bin/bash
#$ -o /storage/essicd/data/HCP/Soroosh/ACAnal/MESim/FPR/logs
#$ -e /storage/essicd/data/HCP/Soroosh/ACAnal/MESim/FPR/logs
#$ -l h_rt=00:30:00
#$ -l h_vmem=4G
#$ -N FPRNP_MEq_t${t_cnt}_rho${rho_cnt}
#$ -r y
#$ -t 1-${nRlz}

cd /home/wmrnaq/ACAnal/Sim/MESim/FPR

.  /etc/profile

#echo t${t_cnt}_rho${rho_cnt}
#rho_list=(0 0 20 50 70 90)
#t_list=(0 100 200 600 1200)
#echo ${SGE_TASK_ID}

#filepath=/storage/essicd/data/HCP/Soroosh/ACAnal/MESim/R_MESim_t${t_cnt}_r${rho_cnt}/
#filename=MESim_t${t_list[${t_cnt}]}_r${rho_list[${rho_cnt}]}_${SGE_TASK_ID}.mat

#echo "${filepath}${filename}"
#if [ -e "${filepath}${filename}" ]
#then
#	echo "exists ${filepath}${filename}"
#	exit
#fi

module add matlab
matlab -nodisplay -nosplash -nodesktop -nojvm -r "rho_cnt=${rho_cnt};t_cnt=${t_cnt};${mat_filename}"
EOF
		qsub $FileName
        done
done

#!/bin/bash
#$ -o /storage/essicd/data/HCP/Soroosh/ACAnal/MESim/FPR/logs
#$ -e /storage/essicd/data/HCP/Soroosh/ACAnal/MESim/FPR/logs
#$ -l h_rt=00:30:00
#$ -l h_vmem=4G
#$ -N FPRNP_MEq_t3_rho1
#$ -r y
#$ -t 1-2000

cd /home/wmrnaq/ACAnal/Sim/MESim/FPR

.  /etc/profile

#echo t3_rho1
#rho_list=(0 0 20 50 70 90)
#t_list=(0 100 200 600 1200)
#echo 

#filepath=/storage/essicd/data/HCP/Soroosh/ACAnal/MESim/R_MESim_t3_r1/
#filename=MESim_t_r_.mat

#echo ""
#if [ -e "" ]
#then
#	echo "exists "
#	exit
#fi

module add matlab
matlab -nodisplay -nosplash -nodesktop -nojvm -r "rho_cnt=1;t_cnt=3;FPR_ME_Sim"

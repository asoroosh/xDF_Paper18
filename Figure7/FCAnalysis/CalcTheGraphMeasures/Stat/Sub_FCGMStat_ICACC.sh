for PP in FPP
do
	for Atlas in CC200 ICA200
	do
FileName="ME_FCGM_${Atlas}_${PP}.sh"
cat > $FileName << EOF
#!/bin/bash
#$ -o /home/wmrnaq/HCP/FC/Mats/Stat/logs
#$ -e /home/wmrnaq/HCP/FC/Mats/Stat/logs
#$ -l h_rt=00:40:00
#$ -l h_vmem=5G
#$ -N ${Atlas}_${PP}_GM_100UR
#$ -r y
#$ -t 1-100

cd ~/HCP/FC/Mats/Stat

.  /etc/profile
module add matlab

matlab -nojvm -nodisplay -nosplash -nodesktop -r "PP='${PP}';Atlas='${Atlas}';FCGMStat;"
EOF
                qsub $FileName
done 
	done

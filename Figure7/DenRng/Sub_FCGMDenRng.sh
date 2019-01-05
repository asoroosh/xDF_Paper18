for PP in FPP
do
	for Atlas in Yeo ICA200 Power Gordon 
	  do
		for DenID in `seq 1 8`
		do
FileName="ME_FCGM_Den${DenID}_${Atlas}_${PP}.sh"
cat > $FileName << EOF
#!/bin/bash
#$ -o /home/wmrnaq/HCP/FC/Mats/DenRng/logs
#$ -e /home/wmrnaq/HCP/FC/Mats/DenRng/logs
#$ -l h_rt=00:40:00
#$ -l h_vmem=4G
#$ -N FCGM_Den${DenID}_${Atlas}_${PP}_100UR
#$ -r y
#$ -t 1-100

cd ~/HCP/FC/Mats/DenRng

.  /etc/profile
module add matlab

matlab -nojvm -nodisplay -nosplash -nodesktop -r "PP='${PP}';DenID=${DenID};Atlas='${Atlas}';FCGMDenRng;"
EOF
                qsub $FileName
		done
	done 
done

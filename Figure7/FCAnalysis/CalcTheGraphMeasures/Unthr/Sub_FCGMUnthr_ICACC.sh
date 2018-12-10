for PP in FPP
do
	for Atlas in ICA200 CC200
	  do
FileName="ME_FCGM_${Atlas}_${PP}.sh"
cat > $FileName << EOF
#!/bin/bash
#$ -o /home/wmrnaq/HCP/FC/Mats/Unthr/logs
#$ -e /home/wmrnaq/HCP/FC/Mats/Unthr/logs
#$ -l h_rt=00:60:00
#$ -l h_vmem=4G
#$ -N ${Atlas}_${PP}_GM_100UR
#$ -r y
#$ -t 1-100

cd ~/HCP/FC/Mats/Unthr

.  /etc/profile
module add matlab

matlab -nojvm -nodisplay -nosplash -nodesktop -r "PP='${PP}';Atlas='${Atlas}';FCGMUnthr;"
EOF
                qsub $FileName
done 
	done

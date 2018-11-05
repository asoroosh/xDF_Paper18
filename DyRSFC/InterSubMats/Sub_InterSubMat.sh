for Atlas in YEO
    do
FileName="InterSubMat_${Atlas}.sh"
cat > $FileName << EOF
#!/bin/bash
#$ -o /home/wmrnaq/DyConn/InterSubMats/logs
#$ -e /home/wmrnaq/DyConn/InterSubMats/logs
#$ -l h_rt=01:30:00
#$ -l h_vmem=5G
#$ -N ${Atlas}_DyConnInterSub_100UR
#$ -r y
#$ -t 1-5000

cd /home/wmrnaq/DyConn/InterSubMats

.  /etc/profile
module add matlab

matlab -nojvm -nodisplay -nosplash -nodesktop -r "Atlas='${Atlas}';InterSubMat"
EOF
                qsub $FileName
done 

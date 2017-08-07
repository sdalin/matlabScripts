#!/bin/sh
# request bourne shell as shell for job
#$ -S /bin/sh
#$ -V
#$ -cwd
#$ -m e
#$ -M sdalin@mit.edu

cd ./
module load jre/1.8.0-77
matlab -nodisplay -r "loop"





#!/bin/bash

STRAND_LIST=~/GitCode/MitoImpute/scripts/INFORMATION_LISTS/b37_platforms.txt

for i in `cat ${STRAND_LIST}`; do
	qsub -v MtPlatforms=${i} ~/GitCode/MitoImpute/scripts/BASH/ThousandGenomes_raijin.sh
	echo "JOB FOR ${i} SUBMITTED"
done

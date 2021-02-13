#$ -S /bin/sh

# gets working directory & fastq_ids
# submits alignment job and variant calling job
# assumes a queue and a parallel environment

BASEDIR=$(dirname $0)

# submits alignment job
qsub -q all.q -N ALN_$2 -pe threaded 12 $BASEDIR/alignment.sh $1 $2  

# submits variant callling 
# waits for alignment job (-hold_jid takes jobid)
qsub -q all.q -N CALL_$2 -pe threaded 3 -hold_jid ALN_$2 $BASEDIR/variantcalling.sh $1 $2 

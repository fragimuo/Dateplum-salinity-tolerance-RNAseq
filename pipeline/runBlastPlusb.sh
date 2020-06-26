#!/bin/bash -l
#SBATCH -p all
#SBATCH -n 8 --mem=150G -t 2-00:00:00
## no mail at the time
###SBATCH --mail-type=ALL

## stop on error
set -e

## be verbose and extend the commands
set -x

## load the modules
module load bioinfo-tools
module load blast

## check the options if any
FMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
EVALUE=1e-5
PROC=8
OPTIONS="-max_target_seqs 20"

## usage
usage(){
echo >&2 \
"
	Usage: $0 [options] <blast command> <fasta file> <index> <out dir>

	Options:
                -f the output format (default to $FMT)
                -e the e-value (default to $EVALUE)
		-p number of threads to use (default $PROC)
		-b blast options (example -task blastn instead of megablast)

        Note:
	     Providing the -f option is broken at the moment
"
	exit 1
}

## get the options
while getopts e:f:p:b: option
do
        case "$option" in
	    f) FMT=$OPTARG;;
	    e) EVALUE=$OPTARG;;
	    p) PROC=$OPTARG;;
	    b) OPTIONS=$OPTARG;;
		\?) ## unknown flag
		usage;;
        esac
done
shift `expr $OPTIND - 1`

## we get two dir and two files as input
if [ $# != 4 ]; then
    echo "This function takes one blast command, one file, one index and one directory as arguments"
    usage
fi

BLAST=
EXT="nhr"
case "$1" in
    blastn)
	;;
    blastp)
	EXT="phr"
	;;
    blastx)
	EXT="phr"
	;;
    tblastn)
	;;
    tblastx)
	;;
    *)
	echo "Unknown blast command. Not one of blastn, blastp, blastx, tblastn, tblastx. Aborting."
	usage
esac
BLAST=$1

if [ ! -f $2 ]; then
    echo "The second argument needs to be the fasta file"
    usage
fi

if [ ! -f $3.$EXT ] && [ ! -f $3.${EXT//hr/al} ]; then
    echo "The third argument needs to be the basename of the BLAST index - including the path"
    usage
fi

if [ ! -d $4 ]; then
    echo "The forth argument needs to be the output directory"
    usage
fi

## run BLAST
echo Blasting
fnam=$(basename $3)_$(basename ${2//.f*a*/.blt})
$1 -db $3 -query $2 -out $4/$fnam -evalue $EVALUE $OPTIONS -num_threads $PROC -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'

##
echo Done



#!/bin/bash

THREADS=1
MINALIGN=40

POSITIONAL=()
while [[ $# -gt 0 ]]
do
	key="$1"
	case $key in
		-a|--annovarfile)
		ANNOVARFILE="$2"
		shift
		shift
		;;
		-b|--bamfiles)
		IFS="," read -ra BAMFILES <<< "$2"
		shift
		shift
		;;
		-t|--threads)
		THREADS="$2"
		shift
		shift
		;;
		-m|--min-alignment-score)
		MINALIGN="$2"
		shift
		shift
		;;
		-h|--help)
			echo "Usage: runScript.sh -a ANNOVAR-file -b COMMA-SEPARATED-LIST-OF-BAM-files -t threads [-m minimum-alignment-score] [-h (help)]"
			exit 0
		shift
		;;
	esac
done
set -- "${POSITIONAL[@]}"

if [ -z "$ANNOVARFILE" ]
then
	echo "Please provide an ANNOVAR-annotated file."
	exit -1
fi

if [ ! -f "$ANNOVARFILE" ]
then
	echo "The ANNOVAR file provided does not exist."
	exit -2
fi

regex='^[0-9]+$'

if ! [[ $THREADS =~ $regex ]]
then
	echo "Please provide an integer value for the number of threads."
	exit -3
fi

if ! [[ $MINALIGN =~ $regex ]]
then
	echo "Please provide an integer value for the minimum alignment score."
	exit -4;
fi

for f in "${BAMFILES[@]}"
do
	if [ ! -f $f ]
	then
		echo "The following BAM file does not exist: $f"
		exit -5
	fi
done

AnnotateBAMStatistics --threads $THREADS --pileup-regions --min-alignment-score $MINALIGN --annovar-file $ANNOVARFILE --bam-files $(printf ",%s" "${BAMFILES[@]}" | sed 's/^,//')

if [[ $? -ne 0 ]]
then
       	echo -e 'Error: An error has been encountered with AnnotateBamStatistics'
       	exit -6
fi
	
exit 0

#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo "Usage: sim_wrapper.sh --out <out_name> --mutrate <float|comma-sep list of floats> --samples <int> --migration <filepath>"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -o|--out) NAME="$2"; shift ;;
        -m|--mutrate) MUTRATE="$2"; shift ;;
        -s|--samples) NUM_SAMPLES="$2"; shift ;;
        -mm|--migration) MIGRATION_MATRIX="$2"; shift ;;

    *) echo "Unknown parameter passed: $1"; echo "Usage: sim_wrapper.sh --out <out_name> --mutrate <float|comma-sep list of floats> --samples <int> --migration <filepath>"; exit 1 ;;
    esac
    shift
done

if [ ! -f "./scripts/simulator.py" ]
then
    echo "Script ./scripts/simulate.py not found. Exiting!"
    exit
fi

if ! command -v reformat.sh &> /dev/null
then
    echo "BBMAP tool reformat.sh could not be found. Exiting!"
    exit
fi

if ! command -v art_illumina &> /dev/null
then
    echo "art_illumina short read simulator could not be found. Exiting!"
    exit
fi

MOUSE="MMUS1469"
REF="BC10v0"
TAG="MG_120419"
ALL="ALL_TISSUES"
FA_ALL="${NAME}"
PREFIXALL="${MOUSE}_${ALL}_${REF}_${TAG}"


# you will likely want to create and activate a conda env using the provided yaml file
# conda env create -f env/simulate.yaml
# conda activate simulate

./scripts/simulator.py ${NAME} ${MUTRATE} ${NUM_SAMPLES} ${MIGRATION_MATRIX}

outputdir="sim_results_${NAME}/"
mkdir ${outputdir}
mv ${NAME}* ${outputdir}

all_dir="${outputdir}all_samples_fastq/"
mkdir ${all_dir}
art_illumina -ss MSv1 -amp -p -na -i ${outputdir}${FA_ALL}.fa -l 175 -f 1000 -o ${all_dir}${PREFIXALL}_R >> ${NAME}.log 2>&1
for file in ${all_dir}*.fq; do mv "$file" "${file%.fq}.fastq"; done

tissues_column=$(cut -f 2 "${outputdir}${NAME}_tissues.tsv" | tail -n +2)
IFS=$'\n' read -d '' -ra TISSUES <<< "$(echo "$tissues_column" | sort | uniq)"
NTISSUES=${#TISSUES[@]}

tissue_dir="${outputdir}tissue_specific_fastq/"
mkdir ${tissue_dir}
for (( i=0; i<${NTISSUES}; i++ )); do
    FA_PREFIX="${NAME}_${TISSUES[$i]}"
    PREFIX="${MOUSE}_${TISSUES[$i]}_${REF}_${TAG}"
    art_illumina -ss MSv1 -amp -p -na -i ${outputdir}${FA_PREFIX}.fa -l 175 -f 1000 -o ${tissue_dir}${PREFIX}_R >> ${NAME}.log 2>&1
done
for file in ${tissue_dir}*.fq; do mv "$file" "${file%.fq}.fastq"; done


### Un-comment below to generate fastq files for each sample from the fasta file with all samples

#mkdir temp_sample_fa/
#awk '/^>/ { OUT=substr($0,2) ".fa" } { print >> "temp_sample_fa/"OUT }' ${outputdir}${FA_ALL}.fa
#samples=$(ls temp_sample_fa/)
#sample_dir="${outputdir}sample_specific_fastq/"
#mkdir ${sample_dir}
#for file in ${samples}; do
#    ID="${file%.fa}"
#    PREFIXSAMPLE="${MOUSE}_${ID}_${REF}_${TAG}"
#    art_illumina -ss HS25 -amp -p -na -i "temp_sample_fa/${file}" -l 150 -f 1000 -o ${sample_dir}${PREFIXSAMPLE}_R >> ${NAME}.log 2>&1
#done
#rm -r temp_sample_fa/
#for file in ${sample_dir}*.fq; do mv "$file" "${file%.fq}.fastq"; done

mv ${NAME}.log ${outputdir}



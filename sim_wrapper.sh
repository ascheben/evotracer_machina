#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo "Usage: sim_wrapper.sh --out <out_name> --mutrate1 <float> --mutrate2 <float> --max-indel-size <int> --samples <int> --migration <file_path>"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -o|--out) NAME="$2"; shift ;;
        -m1|--mutrate1) MUTRATE1="$2"; shift ;;
        -m2|--mutrate2) MUTRATE2="$2"; shift ;;
        -i|--max-indel-size) MAX_INDEL_SIZE="$2"; shift ;;
        -s|--samples) NUM_SAMPLES="$2"; shift ;;
        -mm|--migration) MIGRATION_MATRIX="$2"; shift ;;

    *) echo "Unknown parameter passed: $1"; echo "Usage: sim_wrapper.sh --out <out_name> --mutrate1 <float> --mutrate2 <float> --max-indel-size <int> --samples <int> --migration <file_path>" ; exit 1 ;;
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
RAW_TISSUES=$(head -n 1 ${MIGRATION_MATRIX})
IFS=', ' read -r -a TISSUES <<< "${RAW_TISSUES}"
NTISSUES=${#TISSUES[@]}
ALL="ALL_TISSUES"
FA_ALL="${NAME}"
PREFIXALL="${MOUSE}_${ALL}_${REF}_${TAG}"


# you will likely want to create and activate a conda env using the provided yaml file
# conda env create -f env/simulate.yaml
# conda activate simulate

./scripts/simulator.py ${NAME} ${MUTRATE1} ${MUTRATE2} ${MAX_INDEL_SIZE} ${NUM_SAMPLES} ${MIGRATION_MATRIX}

outputdir="sim_results_${NAME}/"
mkdir ${outputdir}
mv ${NAME}* ${outputdir}

# Previous code for downsampling, but will need to be re-written for new variable framework without PREFIX1, etc.
#reformat.sh in=${PREFIX1}.fa out=${PREFIX2}.fa samplerate=0.5 >> ${NAME}.log 2>&1
#reformat.sh in=${PREFIX1}.fa out=${PREFIX3}.fa samplerate=0.1 >> ${NAME}.log 2>&1

all_dir="${outputdir}all_samples_fastq/"
mkdir ${all_dir}
art_illumina -ss HS25 -amp -p -na -i ${outputdir}${FA_ALL}.fa -l 150 -f 1000 -o ${all_dir}${PREFIXALL}_R >> ${NAME}.log 2>&1
for file in ${all_dir}*.fq; do mv "$file" "${file%.fq}.fastq"; done

tissue_dir="${outputdir}tissue_specific_fastq/"
mkdir ${tissue_dir}
for (( i=1; i<${NTISSUES}; i++ )); do
    FA_PREFIX="${NAME}_${TISSUES[$i]}"
    PREFIX="${MOUSE}_${TISSUES[$i]}_${REF}_${TAG}"
    art_illumina -ss HS25 -amp -p -na -i ${outputdir}${FA_PREFIX}.fa -l 150 -f 1000 -o ${tissue_dir}${PREFIX}_R >> ${NAME}.log 2>&1
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



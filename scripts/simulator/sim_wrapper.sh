#!/bin/bash


if [[ $# -eq 0 ]] ; then
    echo "Usage: sim_wrapper.sh --out <out_name> --mutrate <float|comma-sep list of floats> --samples <int> [--migration <filepath>]"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -o|--out) NAME="$2"; shift ;;
        -m|--mutrate) MUTRATE="$2"; shift ;;
        -s|--samples) NUM_SAMPLES="$2"; shift ;;
        -mm|--migration) MIGRATION_MATRIX="$2"; shift ;;

    *) echo "Unknown parameter passed: $1"; echo "Usage: sim_wrapper.sh --out <out_name> --mutrate <float|comma-sep list of floats> --samples <int> [--migration <filepath>]"; exit 1 ;;
    esac
    shift
done

if [ ! -f "./scripts/simulator/simulator.py" ]
then
    echo "Script ./scripts/simulator/simulator.py not found. Exiting!"
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

# you will likely want to create and activate a conda env using the provided yaml file
# conda env create -f env/simulate.yaml
# conda activate simulate

if [[ "$MIGRATION_MATRIX" != "NA" ]]; then

    MOUSE="MMUS1469"
    REF="BC10v0"
    TAG="MG_120419"
    ALL="ALL_TISSUES"
    FA_ALL="${NAME}"
    PREFIXALL="${MOUSE}_${ALL}_${REF}_${TAG}"

    ./scripts/simulator/simulator.py ${NAME} ${MUTRATE} ${NUM_SAMPLES} ${MIGRATION_MATRIX}

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

else

    MOUSE="MMUS1469"
    TISSUE1="tissue1"
    TISSUE2="tissue2"
    TISSUE3="tissue3"
    REF="BC10v0"
    TAG="MG_120419"
    PREFIX1="${MOUSE}_${TISSUE1}_${REF}_${TAG}_${NAME}"
    PREFIX2="${MOUSE}_${TISSUE2}_${REF}_${TAG}_${NAME}"
    PREFIX3="${MOUSE}_${TISSUE3}_${REF}_${TAG}_${NAME}"

    ./scripts/simulator/simulator.py ${NAME} ${MUTRATE} ${NUM_SAMPLES} "NA"

    outputdir="sim_results_${NAME}/"
    mkdir ${outputdir}
    mv ${NAME}_* ${outputdir}
    mv ${NAME}.* ${outputdir}

    mv ${outputdir}/${NAME}.fa ${outputdir}/${PREFIX1}.fa 

    reformat.sh in=${outputdir}${PREFIX1}.fa out=${outputdir}${PREFIX2}.fa samplerate=0.3 >> ${NAME}.log 2>&1
    reformat.sh in=${outputdir}${PREFIX1}.fa out=${outputdir}${PREFIX3}.fa samplerate=0.1 >> ${NAME}.log 2>&1

    all_dir="${outputdir}all_samples_fastq/"
    mkdir ${all_dir}

    art_illumina -ss MSv1 -amp -p -na -i ${outputdir}${PREFIX1}.fa -l 175 -f 1000 -o ${all_dir}${PREFIX1}_R >> ${NAME}.log 2>&1
    art_illumina -ss MSv1 -amp -p -na -i ${outputdir}${PREFIX2}.fa -l 175 -f 1000 -o ${all_dir}${PREFIX2}_R >> ${NAME}.log 2>&1
    art_illumina -ss MSv1 -amp -p -na -i ${outputdir}${PREFIX3}.fa -l 175 -f 1000 -o ${all_dir}${PREFIX3}_R >> ${NAME}.log 2>&1
    for file in ${all_dir}*.fq; do mv "$file" "${file%.fq}.fastq"; done

fi

mv ${NAME}.log ${outputdir}



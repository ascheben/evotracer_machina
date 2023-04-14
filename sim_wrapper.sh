#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo "Usage: sim_wrapper.sh --out <out_name> --mutrate1 <float> --mutrate2 <float> --max-indel-size <int> --samples <int>"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -o|--out) NAME="$2"; shift ;;
        -m1|--mutrate1) MUTRATE1="$2"; shift ;;
        -m2|--mutrate2) MUTRATE2="$2"; shift ;;
        -i|--max-indel-size) MAX_INDEL_SIZE="$2"; shift ;;
        -s|--samples) NUM_SAMPLES="$2"; shift ;;

    *) echo "Unknown parameter passed: $1"; echo "Usage: sim_wrapper.sh --out <out_name> --mutrate1 <float> --mutrate2 <float> --max-indel-size <int> --samples <int>" ; exit 1 ;;
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
TISSUES=("PRL" "HMR" "LGR")
NTISSUES=${#TISSUES[@]}
ALL="ALL_TISSUES"
FA_ALL="${NAME}"
PREFIXALL="${MOUSE}_${ALL}_${REF}_${TAG}_${NAME}"

# example usage: sim_wrapper.sh sim4 0.1 0.1
# you will likely want to create and activate a conda env using the provided yaml file
# conda env create -f env/simulate.yaml
# conda activate simulate

./scripts/simulator.py ${NAME} ${MUTRATE1} ${MUTRATE2} ${MAX_INDEL_SIZE} ${NUM_SAMPLES}
outputdir="sim_results_${NAME}/"
mkdir ${outputdir}
mv ${NAME}*  ${outputdir}

# Previous code for downsampling, but will need to be re-written for new variable framework without PREFIX1, etc.
#reformat.sh in=${PREFIX1}.fa out=${PREFIX2}.fa samplerate=0.5 >> ${NAME}.log 2>&1
#reformat.sh in=${PREFIX1}.fa out=${PREFIX3}.fa samplerate=0.1 >> ${NAME}.log 2>&1

art_illumina -ss HS25 -amp -p -na -i ${outputdir}${FA_ALL}.fa -l 150 -f 1000 -o ${outputdir}${PREFIXALL}_R >> ${NAME}.log 2>&1

for (( i=0; i<${NTISSUES}; i++ )); do
    FA_PREFIX="${NAME}_${TISSUES[$i]}"
    PREFIX="${MOUSE}_${TISSUES[$i]}_${REF}_${TAG}_${NAME}"
    art_illumina -ss HS25 -amp -p -na -i ${outputdir}${FA_PREFIX}.fa -l 150 -f 1000 -o ${outputdir}${PREFIX}_R >> ${NAME}.log 2>&1
done

mv ${NAME}.log ${outputdir}


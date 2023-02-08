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
TISSUE1="tissue1"
TISSUE2="tissue2"
TISSUE3="tissue3"
REF="BC10v0"
TAG="MG_120419"
PREFIX1="${MOUSE}_${TISSUE1}_${REF}_${TAG}_${NAME}"
PREFIX2="${MOUSE}_${TISSUE2}_${REF}_${TAG}_${NAME}"
PREFIX3="${MOUSE}_${TISSUE3}_${REF}_${TAG}_${NAME}"

# example usage: sim_wrapper.sh sim4 0.1 0.1
# you will likely want to create and activate a conda env using the provided yaml file
# conda env create -f env/simulate.yaml
# conda activate simulate

./scripts/simulator.py ${PREFIX1} ${MUTRATE1} ${MUTRATE2} ${MAX_INDEL_SIZE} ${NUM_SAMPLES}
reformat.sh in=${PREFIX1}.fa out=${PREFIX2}.fa samplerate=0.5 >> ${NAME}.log 2>&1
reformat.sh in=${PREFIX1}.fa out=${PREFIX3}.fa samplerate=0.1 >> ${NAME}.log 2>&1

art_illumina -ss HS25 -amp -p -na -i ${PREFIX1}.fa -l 150 -f 1000 -o ${PREFIX1}_R >> ${NAME}.log 2>&1
art_illumina -ss HS25 -amp -p -na -i ${PREFIX2}.fa -l 150 -f 2000 -o ${PREFIX2}_R >> ${NAME}.log 2>&1
art_illumina -ss HS25 -amp -p -na -i ${PREFIX3}.fa -l 150 -f 500 -o ${PREFIX3}_R >> ${NAME}.log 2>&1


NAME="$1"
MUTRATE1="$2"
MUTRATE2="$3"
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

./scripts/simulator.py ${PREFIX1} ${MUTRATE1} ${MUTRATE2}
reformat.sh in=${PREFIX1}.fa out=${PREFIX2}.fa samplerate=0.5
reformat.sh in=${PREFIX1}.fa out=${PREFIX3}.fa samplerate=0.1

art_illumina -ss HS25 -amp -p -na -i ${PREFIX1}.fa -l 150 -f 1000 -o ${PREFIX1}_R
art_illumina -ss HS25 -amp -p -na -i ${PREFIX2}.fa -l 150 -f 2000 -o ${PREFIX2}_R
art_illumina -ss HS25 -amp -p -na -i ${PREFIX3}.fa -l 150 -f 500 -o ${PREFIX3}_R


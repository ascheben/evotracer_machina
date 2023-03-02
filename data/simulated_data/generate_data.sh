# commands used to generate simulated data
# note that file suffix "fq" must be renamed to "fastq" for EvoTraceR
./sim_wrapper.sh --out simsmall --mutrate1 0.1 --mutrate2 0.01 --max-indel-size 8 --samples 20
./sim_wrapper.sh --out simmid --mutrate1 0.1 --mutrate2 0.05 --max-indel-size 8 --samples 50

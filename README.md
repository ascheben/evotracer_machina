# CRISPR Cas9 lineage tracing data analysis with EvoTraceR/MACHINA and simulation framework for CRISPR Cas9 lineage tracing data

This repo accomplishes 2 main tasks:

1. Capable of analyzing experimental data from a cancer Cas9 lineage tracing system with EvoTraceR to then infer and analyze tumor migration histories with MACHINA.

2. Provides a thorough simulation framework that simulates CRISPR Cas9 lineage tracing barcode sequences as FASTQ files to evaluate the efficiency of the EvoTraceR pipeline in collapsing unique barcode inputs into amplicon sequence variants (ASVs) output. The repo also contains optional capabilities to compare EvoTraceR to other commonly used tools for calling mutations from a single CRISPR Cas9 target site, ampliCan and CRISPResso2. It is also possible to generate simulations based on specifying migration probabilities to evaluate the EvoTraceR and MACHINA migration history inference.


## Required software environments for Simulation, EvoTraceR, MACHINA, ampliCan, and CRISPResso2

### Simulation software
The required software to run the Simulation can be installed with conda:
```
conda env create -f env/simulate.yaml
```
Although Cassiopeia should be installed via pip within the conda environment, in my case I needed to install Cassiopeia separately using pip as described in the [docs](https://cassiopeia-lineage.readthedocs.io/en/latest/installation.html).

### EvoTraceR software
The required software to run EvoTraceR can be installed with conda: 
```
conda env create -f env/r_env.yaml
```
There are Trimmomatic and FLASH dependency filepaths that need to be manually updated in `sim_full_pipeline.sh`. It also requires that the EvoTraceR repo is installed (https://github.com/Nowak-Lab/EvoTraceR) in the parent directory and the path to it is correct in `sim_full_pipeline.sh`.

### MACHINA software
The required software to run MACHINA can be installed with conda:
```
conda env create -f env/machina.yaml
```
Due to inconsistencies in operating systems, the environment file method can fail. In this case, the environment can be create manually as follows:
```
conda create -c conda-forge -c bioconda -c etetoolkit -c anaconda -n machina biopython machina ete3 scipy numpy pandas parafly
```
Note that the proprietary Gurobi solver needs a [license](https://www.gurobi.com/academia/academic-program-and-licenses/). An academic license can be acquired from the link and activated via command-line using `grbgetkey`.

### ampliCan software
The required software to run the ampliCan can be installed with conda:
```
conda env create -n amplican_env -f env/amplican_env.yaml
```

### CRISPResso2 software
The required software to run the CRISPResso2 can be installed with conda:
```
conda env create -f env/crispresso2_env.yaml 
```


## Task 1: Migration history analysis of experimental data with EvoTraceR and MACHINA

To analyze experimental data, we need to run EvoTraceR and then MACHINA.

### EvoTraceR

EvoTraceR is dependent on a certain naming scheme that identifies the tissue labels for where the samples were collected. Here is an example of paired-end reads in FASTQ files labeled for 3 tissues which can also be found in data/fastq/:
```
MMUS1469_HMR_BC10v0_MG_120419_R1.fastq
MMUS1469_HMR_BC10v0_MG_120419_R2.fastq
MMUS1469_LGR_BC10v0_MG_120419_R1.fastq
MMUS1469_LGR_BC10v0_MG_120419_R2.fastq
MMUS1469_PRL_BC10v0_MG_120419_R1.fastq
MMUS1469_PRL_BC10v0_MG_120419_R2.fastq
```
The tissue label must be retained in the same position as HMR, LGR, or PRL as above and the R1 or R2 labels specify the paired reads.

We can now run EvoTraceR with:

```
conda activate r_env
./scripts/evotracer/run_evotracer_parallelize.sh <input_dir> <output_dir> <trimmomatic_path> <flash_path>
```

The key EvoTraceR output files are:
```
asv_stat.csv
tree_all_clones.newick
```

### MACHINA

We can run MACHINA with:

```
conda activate machina
./scripts/machina/run_pipeline.sh --infile data/asv_stat.csv --tree data/tree_all_clones.newick --scripts ./scripts/machina/ --prefix myprefix --primary-tissue PRL
```

Only a single tissue can be designated as the primary tumor source tissue. The cutoff value determines how many labels (unique sequence variant and tissue combinations) a clonal populations has to have to be analysed using the faster "OLD" machina algorithm. If no cutoff value is supplied, all populations will be analysed with the "OLD" algorithm.

The key MACHINA output files are:

```
myprefix_cp_output/myprefix_all_results.txt
myprefix_cp_output/myprefix_migration.txt
myprefix_cp_output/myprefix_seeding_topology.txt
myprefix_cp_output/myprefix_selection_original_expansion.txt
myprefix_cp_output/myprefix_selection_original_test.txt
```

### Boostrap analysis

Cassiopeia needs to be installed to execute a bootstrap analysis by resampling the indel matrix. An indel matrix has been provided here in the data directory. We can resample from this to generate a file with a single newick bootstrap tree per line.

```
./scripts/bootstrap.py > data/MMUS1469_100bootstraps.nwk
cd data
split -l 1 MMUS1469_100bootstraps.nwk -d tree_
rename _900 _9 tree_*
rename tree_0 tree_ tree_*
```

Now we can create a series of commands to run the MACHINA pipeline on each of the bootstrap trees:

`seq 0 99| while read l; do echo "./scripts/machina/run_pipeline.sh --infile data/asv_stat.csv --tree data/tree_${l} --scripts ./scripts/machina/ --prefix ${l}" > run.cmd;done`


## Task 2: Simulating Cas9 barcodes to evaluate EvoTraceR/MACHINA or EvoTraceR, ampliCan, and CRISPResso2.

We developed a simulation framework of ground truth amplicon sequences that is compatible with analysis by EvoTraceR, ampliCan, and CRISPResso2 or by EvoTraceR/MACHINA. The `./scripts/simulator.py` script calls [Cassiopeia](https://cassiopeia-lineage.readthedocs.io/en/latest/index.html) and [art](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) to carry out the simulation. The `./scripts/simulator.py` is wrapped for formatting with `./scripts/simulator/sim_wrapper.sh` and the whole simulation is governed by and best used through `./sim_full_pipeline.sh` which has the following input requirements:

Required command line inputs:
```
[--out <out_name>] is the desired output directory name for the simulation. It should not be an already existing directory.
[--mutrate <float|comma-sep list of floats>] is a 10 float comma seperated list of mutation rates for each of the 10 target sites in the CRISPR Cas9 barcode. Mutation rates must be between 0 and 1.
[--samples <int>] is the number of sample barcodes to be simulated and must be less than 10,000 to be compatible with the current Cassiopeia tree downsampling setup.
```

Optional command line inputs:
```
[--evotracer] is an optional flag specifying to run EvoTraceR.
[--amplican] is an optional flag specifying to run ampliCan.
[--crispresso2] is an optional flag specifying to run CRISPResso2.
[--machina] is an optional flag specifying to run MACHINA. If this optional flag is used, then EvoTraceR is automatically run in combination and the other optional flags are invalid.
[--migration <migration_matrix_filepath>] is optional if not running MACHINA, but required if running MACHINA. It specifies the filepath of the csv matrix with tissue labels and probabilities of migration between tissues for the simulation to output tissue labeled data compatible with MACHINA.
```

### Running the full simulation pipeline with optional EvoTraceR, ampliCan, and CRISPResso2 analysis or EvoTraceR/MACHINA analysis

Here are some example inputs dending on which packages in the simulation pipeline are being run:

The simulation can be run alone, in combination with only EvoTraceR only, or in combination with EvoTraceR/MACHINA without limitations imposed on the user specified mutation rate array for the 10 target sites:

To only run the Simulation:
```
./sim_full_pipeline.sh --out simmid --mutrate 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1 --samples 100
```
or
```
./sim_full_pipeline.sh --out simmid --mutrate 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1 --samples 100 --migration data/migration_matrices/true_migration_prob_matrix.csv
```

To run the Simulation and EvoTraceR:
```
./sim_full_pipeline.sh --evotracer --out simmid --mutrate 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1 --samples 100
```

To run the Simulation and EvoTraceR/MACHINA:
```
./sim_full_pipeline.sh --machina --out simmid --mutrate 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1 --samples 100 --migration data/migration_matrices/true_migration_prob_matrix.csv
```

There are three optional flags to run EvoTraceR, ampliCan, and CRISPResso2 following the Simulation. Each flag can be used independently of whether the other two are included, but here is example input using all three analysis packages by specifying all three optional flags (this may take between 20-30 minutes to run due to separate aligment steps for each package):
```
./sim_full_pipeline.sh --evotracer --amplican --crispresso2 --out simmid --mutrate 0,0,0,0,0.1,0,0,0,0,0 --samples 100
```
*Note that ampliCan analysis requires a barcode with only site 5 active due to the ampliCan limitation to one site and the amplicon setup for cutsite 5. CRISPResso2 is also setup only for site 5 for comparison, but it can accomodate more sites in a future design that does not specify the sgRNA in the command line input.

### Parallelization of the full simulation pipeline

There is a `parallel_sim.sh` wrapper script that allows for manual specification of the inputs to `sim_full_pipeline.sh` to then run the pipeline in parallel batches across nodes using ParaFly. This script is best used after having a general understanding of `sim_full_pipeline.sh` due to the need for extensive manual changes of the `sim_full_pipeline.sh` input within `parallel_sim.sh`. The `parallel_sim.sh`` script can be run with:
```
./parallel_sim.sh
```
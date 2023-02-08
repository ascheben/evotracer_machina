# Migration history analysis

This pipeline uses EvoTracR data based on a cancer Cas9 lineage tracing system to infer and analysie tumor
migration histories.

## Running the pipeline on EvoTracR outputs

Install the dependencies:

```
conda env create -f env/machina.yaml
conda activate machina
```

Note that the proprietary Gurobi solver needs a [license](https://www.gurobi.com/academia/academic-program-and-licenses/). An academic license can be acquired from the link and activated via command-line using `grbgetkey`.

Now we can run the pipeline with MACHINA.

`./run_pipeline.sh --infile data/asv_stat.csv --tree data/tree_all_clones.newick --scripts ./scripts/ --prefix myprefix --primary-tissue prostate --cutoff 50`

Only a single tissue can be designated as the primary tumor source tissue. The cutoff value determines how many labels (unique sequence variant and tissue combinations) a clonal populations has to have to be analysed using the faster "OLD" machina algorithm. 

## Outputs
The key outputs are the files shown below.

```
myprefix_cp_output/myprefix_all_results.txt
myprefix_cp_output/myprefix_migration.txt
myprefix_cp_output/myprefix_seeding_topology.txt
myprefix_cp_output/myprefix_selection_original_expansion.txt
myprefix_cp_output/myprefix_selection_original_test.txt
```

## Bootstrap analysis
Cassiopeia needs to be installed ([see instructions](https://cassiopeia.readthedocs.io/en/latest/setup.html)) to execute a bootstrap analysis by resampling the indel matrix. An indel matrix has been provided here in the `data` directory. We can resample from this to generate a file with a single newick bootstrap tree per line.

```
./scripts/bootstrap.py > data/MMUS1469_100bootstraps.nwk
cd data
split -l 1 MMUS1469_100bootstraps.nwk -d tree_
rename _900 _9 tree_*
rename tree_0 tree_ tree_*
```
Now we can create a series of commands to run the MACHINA pipeline on each of the bootstrap trees.
`seq 0 99| while read l; do echo "./run_pipeline.sh --infile data/asv_stat.csv --tree data/tree_${l} --scripts ./scripts/ --prefix ${l}" > run.cmd;done`

# Simulating test data for EvoTraceR
To evaluate the ability of EvoTraceR to detect mutations in the barcode sequence, we can simulated ground truth amplicon sequences. The `./scripts/simulator.py` script calls [Cassiopeia](https://cassiopeia-lineage.readthedocs.io/en/latest/index.html) and [art](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) to carry out the simulation. A wrapper script is available to generate output suitable for EvoTracer.

```
sim_wrapper.sh --out <out_name> --mutrate1 <float> --mutrate2 <float> --max-indel-size <int> --samples <int>
```

The required software to run the wrapper can be installed with conda:

```
conda env create -f env/simulate.yaml
conda activate simulate
```

Although Cassiopeia should be installed via pip within the conda environment, in my case I needed to install Cassiopeia separately using pip as described in the [docs](https://cassiopeia-lineage.readthedocs.io/en/latest/installation.html).


The provided sample simulated data were generated as shown below. Note that these commands lead to non-determistic mutation profiles.

```
./sim_wrapper.sh --out simsmall --mutrate1 0.1 --mutrate2 0.01 --max-indel-size 3 --samples 20
./sim_wrapper.sh --out simmid --mutrate1 0.1 --mutrate2 0.05 --max-indel-size 5 --samples 50
```

The generated FASTQ files can be used as EvoTraceR input. The `asv_stat.csv` is the only output required for evaluation. The evaluation for EvoTraceR are:

* What percentage of unique mutation combinations in a single sequence were identified as independent Amplicon Sequence Variants (ASVs)?
* What percentage of mutations were identified?
* What percentage of mutations had accurate positions and indel characters?

To compare the `asv_stat.csv` file to the ground truth, the simulated `.fa` files or the separate indel matrix, mutations, cut positions files can be used. Tissue1 data always contains all mutations, whereas other tissues contain random subsets of the full FASTQ file for tissue1. Note that gap character for deletions are stripped from the `.fa` files before being output by `simulator.py`. To get `.fa` files with gaps, this final step can be easily modified to additionally produce files with the gaps. Alternatively the [muscle aligner](https://anaconda.org/bioconda/muscle) can be used to align sequences in the `.fa` files and introduce gaps.



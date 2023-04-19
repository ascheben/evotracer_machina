# Migration history analysis

This pipeline uses EvoTracR data based on a cancer Cas9 lineage tracing system to infer and analysie tumor
migration histories.

## Running the pipeline on EvoTracR outputs

Install the dependencies:

```
conda env create -f env/machina.yaml
conda activate machina
```

Due to inconsistencies in operating systems, the environment file method can fail. In this case, the environment can be create manually as follows:

```
conda create -c conda-forge -c bioconda -c etetoolkit -c anaconda -n machina biopython machina ete3 scipy numpy pandas parafly
```

Note that the proprietary Gurobi solver needs a [license](https://www.gurobi.com/academia/academic-program-and-licenses/). An academic license can be acquired from the link and activated via command-line using `grbgetkey`.

Now we can run the pipeline with MACHINA.

`./run_pipeline.sh --infile data/asv_stat.csv --tree data/tree_all_clones.newick --scripts ./scripts/ --prefix myprefix --primary-tissue prostate`

Only a single tissue can be designated as the primary tumor source tissue. The cutoff value determines how many labels (unique sequence variant and tissue combinations) a clonal populations has to have to be analysed using the faster "OLD" machina algorithm. If no cutoff value is supplied, all populations will be analysed with the "OLD" algorithm.

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
To evaluate the ability of EvoTraceR to detect mutations in the barcode sequence, we simulated ground truth amplicon sequences. The `./scripts/simulator.py` script calls [Cassiopeia](https://cassiopeia-lineage.readthedocs.io/en/latest/index.html) and [art](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) to carry out the simulation. A wrapper script is available to generate output suitable for EvoTracer.

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

The script used for evaluation is 'test_evotracer.py'. Set input parameters in 'test_config.py' such as the file paths for the indel character matrix, the cut positions, the mutations, and the asv_stat.csv file output from EvoTraceR. 

## Test EvoTraceR Usage
```
Usage: python test_evotracer.py [options]

Options:
        -h, --help
                Show this help message and exit
        -d, --output-dir
                Directory for saving selected output file results [default "test_output"]
        -a, --output-all-mutations
                Output all simulated and evotracer ASV cut-site mutation mappings to files. Simulated data is saved to "simulated_cut_site_mutation_mappings.txt" and EvoTraceR data is saved to "evotracer_cut_site_mutation_mappings.txt". Each dictionary represents an ASV's set of mutations mapping its cut site to the indel sequence
        -s, --output-unidentified-sim-mutations
                Output all simulated ASV mutations that do not have an exact match or similar match generated by EvoTraceR to the file "unidentified_simulated_mutations.txt"
        -i, --output-matched-sim-mutations
                Output all simulated ASV mutations that have an exact match generated by EvoTraceR to the file "matched_simulated_mutations.txt"
        -m, --map-similar-evo-sim-mutations
                EvoTraceR ASV mutations without an exact match in the simulated data are mapped to the most similar ASV mutations in the simulated data. Similar mutations' positions are shifted by at most the length of the mutated sequence and characters are rotated by at most one character. Results are output to the file "similar_evotracer_simulated_mutation_mappings.txt"
```

To compare the `asv_stat.csv` file to the ground truth, the separate indel matrix, mutations, cut positions files are used. Tissue1 data always contains all mutations, whereas other tissues contain random subsets of the full FASTQ file for tissue1. Note that gap character for deletions are stripped from the `.fa` files before being output by `simulator.py`. To get `.fa` files with gaps, this final step can be easily modified to additionally produce files with the gaps. Alternatively the [muscle aligner](https://anaconda.org/bioconda/muscle) can be used to align sequences in the `.fa` files and introduce gaps.

## Full simulation pipeline with EvoTraceR and Machina inference

With the `simulate`, `r_env`, and `machina` conda environments already installed and validated on the independent wrappers, the entire simulation -> evotracer -> machina pipeline can be run with the `sim_full_pipeline.sh` script. There are some filepaths for installed dependencies that need to be manually updated for running machina in the last section of the script. Note: this pipeline also requires the EvoTraceR-parallelize repo to be installed and initialized in R.

Here is example input to run the script:
```
./sim_full_pipeline.sh --out simmid --mutrate 0.05 --max-indel-size 5 --samples 50 --migration data/true_migration_prob_matrix.csv
```
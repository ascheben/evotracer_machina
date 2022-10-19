# Migration history analysis (EvoTracR with prostate cancer Cas9 lineage tracing system)

## Running the pipeline on EvoTracR outputs

Install the dependencies:

```
conda env create -f env/machina.yml
conda activate machina
```

Note that the proprietary Gurobi solver needs a [license](https://www.gurobi.com/academia/academic-program-and-licenses/). An academic license can be acquired from the link and activated via command-line using `grbgetkey`.

Now we can run the pipeline with MACHINA.

`./run_pipeline.sh --infile data/asv_stat.csv --tree data/tree_all_clones.newick --scripts ./scripts/ --prefix myprefix`

## Outputs




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


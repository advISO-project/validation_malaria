# Running the simulation pipeline
This file shows the commands that need to be run in order to create simulated reads with the config files created in the ```prepare_simulation_run/scripts/prepare_files_for_simulation.ipynb``` Jupyter notebook. For full details about the files and how they were created, please refer to the notebook itself.  

## Step 1: download the simulation pipeline
Git clone the read simulation pipeline to a local path of your choice:  
```
git clone https://github.com/ISOinaBox/pop_var_sim.git
```
For more details on the pipeline and how to run it, please refer to its repository on [GitHub](https://github.com/advISO-project/pop_var_sim).

## Step 2: run the pipeline
THe following command runs the simulation pipeline for the "WT and crt mutant" dataset. The run creates a folder "results", inside which you will find a folder "fastqs" which contains gzipped FASTQ files for simulated paired-end runs. Each file is named according to the data provided in the sample design file. The name indicates the genotype (wt or crt mutant), the primer panel (spec, grc1 or grc2) and the read count, referring to the quantiles of read counts for the three primer panels from real-world data (q25, q50, q75).  
This command uses the "singularity" profile, but will also work with profile "docker".  

```
nextflow run YOUR-INSTALLATION-PATH/pop_var_sim/main.nf \
    --haplotype_manifest simulation_run_input_files/hap_manifest_example.csv \
    --sample_design_file simulation_run_input_files/sample_manifest_example.json\
    --run_name "wt_and_crt_mut" \
    -profile singularity \
```

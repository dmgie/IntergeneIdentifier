# How to use it

To use this program, simply clone the entire repository , and place your reads within the same folder. Or alternatively, place the contents of your folder onto the parent folder of your reads i.e 

```
/
|script1.sh
|workflowt_script.sh
|----|reads
```

The general steps to run a similar pipeline is:

1. Aligning reads to a reference annotated genome file (hisat, kallisto etc). Convert them to bam if the output was not bam before, sort and index them
2. Run the intergene-finder to create a reference+intergenic file as well as a FASTA file containing all the sequences for the intergenic regions
3. Create the various BED files for the regions of tRNA, CDS and intergenic (each of them should have their gene name/locus tag attached)
4. Run samtools depth using those BED files on each read (i.e 3x for all reads, as there are 3 separate BED files/region types we are looking at)
5. Use the depth-add-name binary to add the names to the depth files for ease of access later on. We now have 'named' depth files containing the read count at each base position (as well as which gene that base position belongs to)
6. For each of these named depth files, we run TPM normalisation (using the python script, `TPM_normalise.py`)
7. We join all the named depth files into one (using the `combined_depths.sh` script)
8. We perform Quantile Normalisation on this combined csv/tsv file. From here we calculate the thresholds to use in the seed-and-extend algorithm
9. Run the seed-and-extend algorithm using the calculated thresholds, as well as the intergenic FASTA created in Step 2. This results in a FASTA file containing sequences that were found to be of significance.
10. We run these sequences through BLAST to check for homologs and significant alignments



# Requirements

- Python3
- Bash
- Samtools
- R


# Notes:
It currently only has been tested on Linux (Ubuntu Server to be specific), so further testing is required. The `depth-add-name` script especially has been compiled for linux, so would require re-compiling for other. The Rust script used to compile the `depth-add-name` binary is in the repo: `dmgie/intergene-finder`
The seed_and_extend_IGR_folding.py script requires the exe files of the ViennaRNA package in a specific folder. If you want to use this script you might have to change the path to the exe file in the code. In the workflow.sh the normal seed_and_extend_IGR.py is called, which does not require the ViennaRNA package.


Additionally, the workflow.sh script has not yet been completed to allow for ease-for-use. Each tool/script in the folder can be used independently at this stage on the data.


### TODO: Get full path of each defined file and folder path

THREADS=3
BASE_FOLDER="./" # Absoloute path of current folder
REFERENCE_FOLDER="$BASE_FOLDER/reference_files"
# TODO: Maybe also use "readlink -m 'path'" to remove extra slashes that might appear

################################# Workflow script for intergenic region annotation #################################

######## Section 1: Alignments ###########

# Requirements: GFF/GTF file, fasta file, reads
# Outputs: SAM files
# Uses: hisat
REFERENCE_FASTA_FILE="$REFERENCE_FOLDER/example.fasta"
HISAT_INDEX_NAME="strept"
HISAT_INDEX_FOLDER="$BASE_FOLDER/hisat_index"
SAM_FOLDER="$BASE_FOLDER/data_bam"
BAM_FOLDER="$BASE_FOLDER/data_sam"

mkdir -p HISAT_INDEX_FOLDER
mkdir -p SAM_FOLDER
mkdir -p BAM_FOLDER


# TODO: Make and goto a folder to write these indexes to
hisat2-build $REFERENCE_FASTA_FILE $HISAT_INDEX_FOLDER/$HISAT_INDEX_NAME

mkdir -p SAM_FOLDER
mkdir -p BAM_FOLDER

for f in *fastq.gz ; do 
	hisat2 --no-spliced-alignment --rna-strandness "R" --threads $THREADS \
		-q -x $HISAT_INDEX_FOLDER/$HISAT_INDEX_NAME -U $f -S $SAM_FOLDER/$f.sam; 
done


######## Section 2: Conversion to BAM and sort  ###########
# Requirements: SAM files
# Outputs: BAM files
# Uses: samtools

# NOTE: ${f%.*} means only the basename, without extension
for f in $SAM_FOLDER/*.sam; do
	samtools sort -O bam -o $BAM_FOLDER/${f%%.*}.bam
done

# Make index of all of them
for f in $BAM_FOLDER/*.bam; do
	samtools index $f
done


######## Section 3: Creation of BED files ###########
# Requirements: GFF File
# Outputs: BED files (3), intergenic.fasta (all intergenic sequences), reference+intergenic.gff (contains all entries of original + intergenic regions)
# Uses: custom script
# Repeat for: Genic, intergenic, tRNA

REFERENCE_GFF_FILE=""
GFF_WITH_INTERGENIC="$BASE_FOLDER/reference+intergenic.gff"
INTERGENE_BED="intergene_regions.bed"
TNRA_BED="trna_regions.bed"
GENIC_BED="genic_regions.bed"
BED_FOLDER="$BASE_FOLDER/bed_files"

# Extracting intergenic regions
# TODO: Run the intergene finder (rust binary)
./intergene-finder $REFERENCE_GFF_FILE $REFERENCE_FASTA_FILE

## BED Creation
# Intergenic with name
cat $GFF_WITH_INTERGENIC | grep intergene-finder | awk -F '\t' '{split($9, a, ";"); $9=a[1]; split($9, b, "="); $9=b[2];  print $1 "\t" $4 "\t" $5 "\t" $9}' > $BED_FOLDER/$INTERGENE_BED

# Intergenic with name -> only above certain value; change min=50 to something else
# cat $GFF_WITH_INTERGENIC | grep intergene-finder | awk -F '\t' -v min=50 '{split($9, a, ";"); $9=a[1]; split($9, b, "="); $9=b[2];  if ($5 - $4 > min) {print $1 "\t" $4 "\t" $5 "\t" $9}}' > intergene_regions.bed

# tRNA with name (either use this or above)
cat $REFERENCE_GFF_FILE | grep -E 'tRNA' | sort -k3,3n | awk -F '\t' '$3 == "gene" {print}' | awk -F '\t' '{split($9,a,";"); $9=a[2]; split($9,b," "); $9=b[2]; print "Chromosome\t" $4 "\t" $5 "\t" $9}' | sed 's/"\(.*\)"/\1/' > $BED_FOLDER/$TRNA_BED


# GENIC BED WITH NAME
cat $REFERENCE_GFF_FILE | grep -E 'gene' | sort -k3,3n | awk -F '\t' '$3 == "gene" {print}' | awk -F '\t' '{split($9,a,";"); $9=a[2]; split($9,b," "); $9=b[2]; print "Chromosome\t" $4 "\t" $5 "\t" $9}' | sed 's/"\(.*\)"/\1/' > $BED_FOLDER/$GENIC_BED


######## Section 4: Creating .depth files (depth on regions) ###########
# Requirements: BED files
# Outputs: Depth files for regions defined in BED
# Uses: samtools (inside a custom bash script) 
# NOTE: Uses the -a option in samtools to get all bases

# Repeat for: Genic, intergenic, tRNA
DEPTH_FOLDER="$BASE_FOLDER/depths"
DEPTH_SCRIPT="$BASE_FOLDER/get_depth.sh"
TRNA_DEPTH_FOLDER="$DEPTH_FOLDER/tRNA_depth"
INTERGENIC_DEPTH_FOLDER="$DEPTH_FOLDER/intergenic_depth"
GENIC_DEPTH_FOLDER="$DEPTH_FOLDER/genic_depth"

mkdir -p $DEPTH_FOLDER
mkdir -p $TRNA_DEPTH_FOLDER
mkdir -p $INTERGENIC_DEPTH_FOLDER
mkdir -p $GENIC_DEPTH_FOLDER

bash $DEPTH_SCRIPT -i $BAM_FOLDER -o $TRNA_DEPTH_FOLDER -b $TRNA_BED
bash $DEPTH_SCRIPT -i $BAM_FOLDER -o $INTERGENIC_DEPTH_FOLDER -b $INTERGENE_BED
bash $DEPTH_SCRIPT -i $BAM_FOLDER -o $GENIC_DEPTH_FOLDER -b $GENIC_BED



######## Section 5: Adding names to depth files (creates .depthn) ###########
# Requirements: BED files, .depth files
# Outputs: Depth files with gene names attached to each row/base position
# Uses: Custom script (rust but a binary)

# Repeat for: Genic, intergenic, tRNA
ADD_NAME_SCRIPT="$BASE_FOLDER/depth-add-name"

# TODO: Add output folder option
./$ADD_NAME_SCRIPT -i $TRAN_DEPTH_FOLDER/* -b $TRNA_BED
./$ADD_NAME_SCRIPT -i $INTERGENIC_DEPTH_FOLDER/* -b $INTERGENE_BED
./$ADD_NAME_SCRIPT -i $GENIC_DEPTH_FOLDER/* -b $GENIC_BED

# TODO: Delete files afterwards

######## Section 6: TPM Normalising the files ###########
# Requirements: .depthn (named) files, BED files, 
# Outputs: TPM normalised .depthn.norm files
# Uses: Python script

# Repeat for: Genic, intergenic, tRNA
TPM_SCRIPT="$BASE_FOLDER/TPM_normalise.py"

# Place in the same folder as normal named depth
python3 TPM_normalise.py -b $TRNA_BED -i $TRNA_DEPTH_FOLDER/* -o $TRNA_DEPTH_FOLDER
python3 TPM_normalise.py -b $INTERGENIC_BED -i $INTERGENIC_DEPTH_FOLDER/* -o $INTERGENIC_DEPTH_FOLDER
python3 TPM_normalise.py -b $GENIC_BED -i $GENIC_DEPTH_FOLDER/* -o $GENIC_DEPTH_FOLDER




######## Section 7: Quantile normalisation and plots ###########
# Requirements: .depthn.norm files (optional: .depthn files to get plots before normalisation - but slower)
# Outputs: PNG of boxplots, CSV files containing all quantile normalised counts
# Uses: R script

# Repeat for: Genic, intergenic, tRNA

# NOTE: tRNA give command line options "n" "y" when running so we get strain-specific medians (as well as seed-and-extend thresholds)
# TOOD: ^ requires us to parse the output of R to get these values -> could save it here in a variable to be used for the python script below
# NOTE: Maybe skip plots / csv stuff? For most of them?

# NOTE: This is too slow, maybe use the column combining script to create a tsv and read it in afterwards into R using fread








######## Section 8: Analysis of quantile normalisde reads ###########
# Requirements: Quantile normalised reads (csv? with gene_names) for INTERGENIC regions, intergenic fasta file
# Outputs: Possible significant intergenic sequences
# Uses: Python script
INTERGENIC_FASTA="$REFERENCE_FOLDER/intergenic.fasta"

# NOTE: Example, change a few things
python3 seed_and_extend_IGR.py -s $INTERGENIC_FASTA -d asdajsodi -c 'r' -t 1860 -w 21 -e 35 -g 300 -o intergenic_found.fasta



######## Section 9: BLAST homology analysis ###########
# Requirements: Databases (rfam, rnacentral), sequences from seed-and-extend
# Outputs: Significant hits (if found)
# Uses: blastn

RNACENTERAL_DB="$BASE_FOLDER/databases/rnacentral_db"
RFAM_DB="$BASE_FOLDER/databases/rfam_db"

# NOTE: blastn-short due to usually short sequences
blastn -task blastn-short -query intergenic_found.fasta -db rnacentral_db
blastn -task blastn-short -query intergenic_found.fasta -db rfam_db




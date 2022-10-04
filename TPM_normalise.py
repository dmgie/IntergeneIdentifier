import argparse

argparser = argparse.ArgumentParser(description='Normalise TPM values')
argparser.add_argument('-i', '--input', help='Input file - can be multiple. If a whole folder is required, do something like folder/*', required=True, nargs='+')
argparser.add_argument('-b', '--bed', help='Bed File', required=True)
argparser.add_argument('-o', '--output', help='Output folder. If not defined will write it to the current folder', default='./')
# argparser.add_argument('-o', '--output', help='Output file', required=True)
args = argparser.parse_args()

# Stuff from BED files to be used on inputs
gene_dict = {}
outputnames = []

# read tab delimited BED file (once)
with open(args.bed, 'r') as f:
    for line in f:
        # Create a dictionary of gene lengths and their name
        columns = line.split('\t')
        start = int(columns[1])
        end = int(columns[2])
        gene = columns[3]
        gene_dict[gene] = end - start



# For each file in args.input (incase multiple files are given)
for file in args.input:
    print(f"Normalising {file} using TPM...")
    # Reset these on each file
    gene_length = 73 / 1000
    desc = []
    rpk = []
    scaling_factor = 1

    # Read the file and get the gene length
    with open(file, "r") as f:
        outputnames.append(file + ".norm")
        # Go through each line
        for line in f:
            # Split line into columns
            columns = line.split("\t")
            chrom = columns[0]
            pos = columns[1]
            count = columns[2]
            # If gene names exist then good, otherwise empty character
            try:
                gene = columns[3]
                gene_length = int(gene_dict.get(gene)) / 1000
            except IndexError:
                print("No gene names found")
                gene = ""
                gene_length = 73 / 1000 # Default to tRNA length
            # Calculate rpk
            rpk.append(float(count) / gene_length)
            desc.append(f"{chrom}\t{pos}\t{gene}")
        scaling_factor = sum(rpk) / 1000000
    
    # Calculate tpm after having calculated scaling factor above
    tpm = [rpk / scaling_factor for rpk in rpk]


    # TODO: Allow defining an output folder
    # Write each file
    print(f"...writing normalised values to {file}.depthn.norm")
    with open(f"{args.output}{file}.norm", "w") as f:
        for i in range(len(tpm)):
            d = desc[i].split("\t")
            f.write(f"{d[0]}\t{d[1]}\t{tpm[i]}\t{d[2]}")






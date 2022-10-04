## Install if neccessary
#if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("preprocessCore")
# library(qsmooth)

library(preprocessCore)
library(data.table)

# THIS SCRIPT REQUIRES THE DEPTH FILES WITH NAMES 
# THIS SCRIPT ALSO REQUIRES TPM NORMALISED FILES IN THE SAME FORMATA

# TODO: Make into a list, and run run_analysis for each pattern in the list 
# Make it a tuple so its (pattern, strainname), so we can give i.e ("c-e", "Strain1")
# strain_specific_patterns <- c("c-e", "f-h")

# TODO: Get overall median of entire matrix -> used for seed-and-extend extension threshold

# TODO: lapply + do.call slow, use better methods i.e rbindrow/cbindrow?

# Pattern of names, what to append to name, and whether to do the before calulcation
run_analysis <- function(pattern, strain, before="n", median="n") {

	# Before TPM normalised files -> according to pattern
	if (before == "y")  {
		temp <- list.files(pattern=paste0("*[", pattern ,"].*.depthn$"), full.names=TRUE)
		beforeTPM <- lapply(temp, fread, select=c(3), sep="\t")
		names(beforeTPM) <- temp
		bTPM <- do.call("cbind", beforeTPM)

		# Boxplot before TPM
		png(file=paste0("./beforeTPM_", strain, ".png"),
		width=1200, 
		height=800)
		boxplot(bTPM)
		dev.off()
	}

	print("Loading files in")
	# This one loads in after TPM normalised data
	temp <- list.files(pattern=paste0("*[", pattern ,"].*.depthn.norm$"), full.names=TRUE)
	afterTPM <- lapply(temp, fread, select=c(3), sep="\t", header=FALSE, showProgress=TRUE)
	# names(afterTPM) <- temp
	aTPM <- do.call("cbind", afterTPM)

	# Get the gene each base position (row) belongs to
	print("Loading gene_names in")
	gene_names <- fread(temp[1], sep="\t", select=c(4))
	# Make column names better
	columnNames <- gsub(temp, pattern=".bam.*", replacement="")


	# Boxplot after TPM (before quantile normalisation)
	print("Creating the first boxplot")
	png(file=paste0("./afterTPM_", strain, ".png"),
	width=1200, 
	height=800)
	boxplot(aTPM)
	dev.off()

	# Boxplot after (after quantile normalisation)
	normalised <- normalize.quantiles(as.matrix(aTPM))
	png(file=paste0("./afterQN_", strain, ".png"),
	width=1200, 
	height=800)
	boxplot(normalised)
	dev.off()

	dimnames(normalised)[1] <- gene_names
	colnames(normalised) <- columnNames
	write.csv(normalised, paste0("./afterQN_", strain, ".csv"))

	# Get median for every gene
	# from second column onwards (since first is the gene_name)
	if (median == "y") {
		per_gene_median <- aggregate(normalised[,2:ncol(normalised)], list(normalised$gene_name), median)
		per_sample_median <- apply(per_gene_median[,2:ncol(per_gene_median)], 2, median)
		per_strain_median <- median(per_sample_median)
		print(per_strain_median)
	}
}

args = commandArgs()
run_analysis("c-e","Strain1", args[1], args[2])
run_analysis("f-h","Strain2", args[1], args[2])

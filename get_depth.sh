# This is to calculate the depth for each base, at each tRNA region.


# outFolder="."
# bedFile=""
# inputFolder="."
while getopts ":h:i:o:b:" option; do
   case $option in
      h) # display Help
         echo "This is run by using ./get_depth.sh in the folder with the bam files. The -o can be used to specify the output folder"
         exit;;
	  i)
		inputFolder=${OPTARG}
		;;
      o) 
		 outFolder=${OPTARG}
		 ;;
      b) 
		 bedFile=${OPTARG}
		 ;;
   esac
done


for f in ${inputFolder}/*.bam; do
	echo "Running depth for file $f" 
	# echo "samtools depth -a -H -b ${bedFile} -o $outFolder$(basename $f).depth $f;"
	samtools depth -a -b ${bedFile} -o $outFolder$(basename $f).depth $f;
	echo "Completed running for $f"
done

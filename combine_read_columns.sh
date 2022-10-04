# This takes all named ".depth" files and makes it look really nice
# Use via ./combined_read_colums.sh *.depth > allCounts.tsv

# Only takes in named files
firstFile=$1
numinputs=$((4 * $# - 1))
# Find a way to make this tab-delimited
# header="basepos $1 gene_name $@"
header="basepos $1 gene_name $@"
echo "$header | sed 's/ /\t/g'"
# 2,3,4 are the basenumber, depth, name and then counts for each
# Then just take the count column from each one
paste -d'\t' "$@"  | cut -f "2,3,4,$(seq -s ',' 7 4 $numinputs)"

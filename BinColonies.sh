#READ: takes in three arguments
#1: Directory name containing the sample identifier
#2: First read set 
#3: Second read set

#Path for MaxBin program
maxbin="/project/biocomplexity/isentry/src/MaxBin-2.2.7/run_MaxBin.pl"

#Contigs directory
contigs="Output/*.contigs.fasta"
#Change prefix as necessary
prefix="$(echo $1 | cut -d"_" -f 1,3)"
f1=$2
f2=$3
mkdir tmp_dir
mkdir MaxBin
$maxbin -thread 8 -contig $contigs -reads $f1 -reads2 $f2 -out tmp_dir/$prefix.bin  
mv tmp_dir/$prefix.bin.*.fasta MaxBin/
rm -r tmp_dir

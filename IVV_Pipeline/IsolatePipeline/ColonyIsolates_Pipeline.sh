#LOAD MODULES: gcc MUST be loaded before diamond
module load gcc
module load diamond

#CHECK REQUIRED PROGRAM PATHS
check_Mash=$(which mash)
check_Spades=$(which spades.py)

#Database paths for Diamond and Patric (Mash)
cardDB="/project/biocomplexity/isentry/ref_data/card/card_protein_variant_DB.dmnd"
vfdbDB="/project/biocomplexity/isentry/ref_data/vfdb/VFDB_protein_DB.dmnd"
patricDB="/project/biocomplexity/isentry/ref_data/mash/patric_all.msh"
patricMapping="/project/biocomplexity/isentry/ref_data/mash/patric_genomes_names.txt"

#Headers for diamond and mash output
dmndHeaders="qseqid sseqid  length  evalue  bitscore    stitle"
mashHeaders="gene-ID    species  distance   p-value shared-hashes"

#Number of threads
threads="8"

#mash threshold: Sets maximum distance value for output a mash pair
mash_threshold="0.05"

#Diamond: evalue threshold
evalue="0.000001"

#Check that database files exists
if [ ! -f $cardDB ] 
then
    echo "card-diamond database does not exist:"$cardDB
    exit 1
fi
if [ ! -f $vfdbDB ]
then
    echo "vfdb-diamond database does not exist:"$vfdbDB
    exit 1
fi
if [ ! -f $patricDB ]
then
    echo "patric-mash database does not exist:"$patricDB
    exit 1
fi
if [ ! -f $patricMapping ]
then
    echo "patric-mapping file does not exist:"$patricMapping
    exit 1
fi

#Create Output directory if it doesn't exist
if [ ! -d $outdir ]
then
    mkdir $outdir
fi

#File and Output Paths
outdir="Output"

#Prefix for the output files
prefix=$outdir"/output"
#Get forward and reverse read filenames
f1=$1
f2=$2

#Run SPAdes to assemble contigs: Assume two files for forward and reverse reads
spades.py --meta -t $threads -1 $f1 -2 $f2 -o tmp_dir 
mv tmp_dir/contigs.fasta $outdir
rm -r tmp_dir

contigs=$outdir"/contigs.fasta"

#Run Diamond: card
echo $dmndHeaders > $prefix".card"
diamond blastx --outfmt 6 $dmndHeaders --db $cardDB --evalue $evalue --query $contigs >> $prefix".card" 

#Run Diamond: vfdb
echo $dmndHeaders > $prefix".vfdb"
diamond blastx --outfmt 6 $dmndHeaders --db $vfdbDB --evalue $evalue --query $contigs >> $prefix".vfdb"

#Run mash: patric
#TODO: Replace this with a function that loops over possible distance values finding the minimum value
#TODO: Create a python script that calls mash as a subprocees and parse through the file using this script
mash dist -d $mash_threshold $patricDB $contigs | cut -f 1,3,4,5  > tmp_mash.txt 
if [ ! -s "tmp_mash.txt" ]
then
    echo "Mash failed at threshold $mash_threshold" > mash_failed.txt
else
    grep -F -w "$(cut -f 1 tmp_mash.txt | cut -d'/' -f 6)" $patricMapping > tmp_patric.txt 
    echo $mashHeaders > $prefix".patric"
    paste tmp_patric.txt tmp_mash.txt | cut -f 1,2,4,5,6 >> $prefix".patric" 
    rm tmp_patric.txt
fi 
rm tmp_mash.txt

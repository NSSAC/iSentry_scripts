#Read first!

#LOAD MODULES: gcc MUST be loaded before diamond
module load gcc
module load diamond

#CHECK REQUIRED PROGRAM PATHS
check_Mash=$(which mash)

#File and Output Paths
#FastqPath="TestFastq/*"
outdir="IVV_BinReports"
reportFile=$outdir"/ContigsReport.txt"

#Database paths for Kraken2, Diamond, and Patric (Mash)
cardDB="/project/biocomplexity/isentry/ref_data/card/card_protein_variant_DB.dmnd"
vfdbDB="/project/biocomplexity/isentry/ref_data/vfdb/VFDB_protein_DB.dmnd"
patricDB="/project/biocomplexity/isentry/ref_data/mash/patric_all.msh"
patricMapping="/project/biocomplexity/isentry/ref_data/mash/patric_genomes_names.txt"

#Headers for diamond and mash output
dmndHeaders="qseqid sseqid  length  evalue  bitscore    stitle"
mashHeaders="gene-ID   species distance   p-value"

#mash threshold: Sets maximum distance value for output a mash pair
#TODO: Change this to a function that can be looped over, 
#checking if the distance results in an empty file. 
#If it does, increase the maximum distance
mash_threshold="0.04"

#Blast evalue
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

bindir="IVV_Report_Output"

#loop through and create a report for each bin made by MaxBin
for bin in $bindir/maxbin.*.fasta;
do
    #Prefix for bin files put into the outdirectory specified above
    bin_prefix=$outdir"/"$(echo $bin | cut -d"/" -f 2 | cut -d"." -f 1,2) 

    #Run diamond with CARD database
    echo $dmndHeaders > $bin_prefix".card"
    diamond blastx --outfmt 6 $dmndHeaders --db $cardDB --evalue $evalue --query $bin >>  $bin_prefix".card"

    #Run diamond with VFDB database
    echo $dmndHeaders > $bin_prefix".vfdb"
    diamond blastx --outfmt 6 $dmndHeaders --db $vfdbDB --evalue $evalue --query $bin >> $bin_prefix".vfdb"

    #Hash the bin file against a sketch of the patric genomoe dataabase
    #TODO: Make this a function with the specifications from above
    echo $mashHeaders > $bin_prefix".patric"
    mash dist -d $mash_threshold $patricDB $bin | cut -f 1,3,4,5  > tmp_mash.txt
    if [ ! -s "tmp_mash.txt" ]
    then
        head tmp_mash.txt
        echo "Mash from $bin at threshold $mash_threshold is empty"
    else    
        grep -F -w "$(cut -f 1 tmp_mash.txt | cut -d'/' -f 6)" $patricMapping > tmp_patric.txt 
        paste tmp_patric.txt tmp_mash.txt | cut -f 1,2,4,5,6 > $bin_prefix".patric"
        rm tmp_patric.txt
    fi
    rm tmp_mash.txt
done

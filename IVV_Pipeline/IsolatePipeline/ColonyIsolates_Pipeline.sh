#LOAD MODULES: gcc MUST be loaded before diamond
module load gcc
module load diamond

#Diamond databases for card and vfdb
cardDB="/project/biocomplexity/isentry/ref_data/card/card_protein_variant_DB.dmnd"
vfdbDB="/project/biocomplexity/isentry/ref_data/vfdb/VFDB_protein_DB.dmnd"

#Number of threads
threads="8"

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

#runspades

#Run Diamond: card and vfdb

#Run Mash

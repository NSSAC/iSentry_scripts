#LOAD MODULES: gcc MUST be loaded before diamond
module load gcc
module load diamond

maxbin="/project/biocomplexity/isentry/src/MaxBin-2.2.7/run_MaxBin.pl"

#add absolute path to run files
abs_path="/scratch/cc8dm/iSentry_scripts/IVV_Pipeline/MetagenomicPipeline/"

#Flags for skipping parts of the pipeline
runSpades="1"
runKraken="1"
runMashScreen="1"
runDiamond="1"
runCheckm="1"

#Database paths for Kraken2, Diamond, and Patric (Mash)
patricDB="/project/biocomplexity/isentry/ref_data/mash/patric_all.msh"
patricMapping="/project/biocomplexity/isentry/ref_data/mash/patric_genomes_names.txt"

#Headers for diamond and mash output
mashHeaders="gene-ID   species distance   p-value"

#mash threshold: Sets maximum distance value for output a mash pair
#TODO: change, probably unnecessary
mash_threshold="0.25"

#Input arguments
threads="1"
sample="none"
outdir="SampleOutput"
forward="none"
reverse="none"
card="0"
vfdb="0"
while getopts ":t:cvs:o:f:r:" arg; do
    case $arg in
        t)
            threads=${OPTARG}
            ;;
        c)
            card="1"
            ;;
        v)
            vfdb="1"
            ;;
        s)
            sample=${OPTARG}
            ;;
        o)
            outdir=${OPTARG}
            ;;
        f)
            forward=${OPTARG}
            ;;
        r)
            reverse=${OPTARG}
            ;;
        \?)
            echo "Invalid option: -$OPTARG" 1>&2
            exit -1
            ;;
    esac
done

#Check that reads file is provided and exists
if [[ $forward == "none" ]]
then
    echo "Provide a reads file: [ -f FILE ]" 1>&2
    exit -1
fi
if [[ ! -f $forward ]]
then
    echo "$forward does not exist" 1>&2
    exit -1
fi

#Create Output directory if it doesn't exist
if [[ ! -d $outdir ]]
then
    mkdir $outdir
fi
cd $outdir

#Run Spades
if [[ "$runSpades" == "1" ]]
then
    sh "$abs_path"RunSpades_Metagenomic.sh -f $forward -r $reverse -t $threads
fi
if [[ ! -f tmp_dir/contigs.fasta ]] && [[ ! -f contigs.fasta  ]]
then
    echo "Spades failed for $forward" 1>&2
    exit -1
elif [[ -f contigs.fasta  ]] 
    echo "Spades has already been run, contigs.fasta exists" 1>&2
then
else
    mv tmp_dir/contigs.fasta .
    rm -r tmp_dir
fi

#Run Kraken on contigs file
if [[ "$runKraken" == "1" ]]
then
    sh "$abs_path"RunKraken.sh -t $threads -b -r -g -s -f contigs.fasta -i $sample   
fi

#Run MashScreen on contigs file
if [[ "$runMashScreen" == "1" ]]
then
    sh "$abs_path"SetupMash_Screen.sh -t $threads -f contigs.fasta -s $sample 
fi

#Run Diamond on contigs file
if [[ "$runDiamond" == "1" ]] 
then
    if [[ "$card" == "1" ]]
    then
        sh "$abs_path"RunDiamond_Metagenomic.sh -c -f contigs.fasta -s $sample
    fi
    if [[ "$vfdb" == "1" ]]
    then
        sh "$abs_path"RunDiamond_Metagenomic.sh -v -f contigs.fasta -s $sample
    fi
fi

#Incorporate Checkm script
if [[ "$runCheckm" == "1" ]] 
then
    sh RunCheckm.sh
fi

#TODO:run MaxBin
mkdir tmp_dir
$maxbin -thread $threads -contig $outdir/contigs.fasta -reads $f1 -reads2 $f2 -out tmp_dir/maxbin
mv tmp_dir/maxbin.*.fasta $outdir
rm tmp_dir/*
rm -r tmp_dir

#TODO:Run kraken on binned files

#TODO: Run MashScreen on binned files

#TODO: Section below is probably redundant, keep correct and remove incorrect parts
#loop through and create a report for each bin made by MaxBin
for bin in $outdir/maxbin.*.fasta;
do
    #Prefix for each bin: Includes path from output directory specified above
    bin_prefix=$outdir"/"$(echo $bin | cut -d"/" -f 2 | cut -d"." -f 1,2) 

    #Run Diamond using CARD database
    echo $dmndHeaders > $bin_prefix."card"
    diamond blastx --outfmt 6 $dmndHeaders --db $cardDB --evalue $evalue --query $bin >> $bin_prefix".card"

    #Run Diamond using VFDB database
    echo $dmndHeaders > $bin_prefix".vfdb"
    diamond blastx --outfmt 6 $dmndHeaders --db $vfdbDB --evalue $evalue --query $bin >> $bin_prefix".vfdb"

    #Mash the fasta files versus a sketch of the entire patric genome
    #TODO: replace this step with a function that loops over possible thresholds
    #TODO: replace this with a python script that calls mash as a subprocess and checks the file
    #TODO: use MashScreen instead of Mash for the metagenomic samples
    #TODO: remove bin processing
    mash dist -d $mash_threshold $patricDB $bin | cut -f 1,3,4,5  > tmp_mash.txt
    if [[ ! -s "tmp_mash.txt" ]]
    then
        echo "Mash from $bin at threshold $mash_threshold is empty"
    else
        grep -F -w "$(cut -f 1 tmp_mash.txt | cut -d'/' -f 6)" $patricMapping > tmp_patric.txt
        echo $mashHeaders > $bin_prefix".patric"
        paste tmp_patric.txt tmp_mash.txt | cut -f 1,2,4,5,6 >> $bin_prefix".patric"
        rm tmp_patric.txt
    fi
    rm tmp_mash.txt
done

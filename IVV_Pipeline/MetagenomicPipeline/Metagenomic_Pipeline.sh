#LOAD MODULES: gcc MUST be loaded before diamond
module load gcc
module load diamond

maxbin="/project/biocomplexity/isentry/src/MaxBin-2.2.7/run_MaxBin.pl"

#Database paths for Kraken2, Diamond, and Patric (Mash)
krakenBDB="/project/biocomplexity/isentry/ref_data/kraken2/bacteria"
krakenVDB="/project/biocomplexity/isentry/ref_data/kraken2/viral"
cardDB="/project/biocomplexity/isentry/ref_data/card/card_protein_variant_DB.dmnd"
vfdbDB="/project/biocomplexity/isentry/ref_data/vfdb/VFDB_protein_DB.dmnd"
patricDB="/project/biocomplexity/isentry/ref_data/mash/patric_all.msh"
patricMapping="/project/biocomplexity/isentry/ref_data/mash/patric_genomes_names.txt"

#Headers for diamond and mash output
dmndHeaders="qseqid sseqid  length  evalue  bitscore    stitle"
mashHeaders="gene-ID   species distance   p-value"

#mash threshold: Sets maximum distance value for output a mash pair
mash_threshold="0.25"

#Diamond evalue threshold
evalue="0.000001"

#Input arguments
threads="1"
sample="none"
outdir="SampleOutput"
forward="none"
reverse="none"
card="0"
vfdb="0"
while getopts "tcvs:of:r" arg; do
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

#Create Output directory if it doesn't exist
if [ ! -d $outdir ]
then
    mkdir $outdir
fi
cd $outdir

sh ../RunSpades_Metagenomic.sh 

#Run kraken on raw fastq input files
kraken2 --db $krakenBDB --threads $threads --paired $f1 $f2 --output $prefix".bacteria.output" --report $prefix".bacteria.report"
kraken2 --db $krakenVDB --threads $threads --paired $f1 $f2 --output $prefix".viral.output" --report $prefix".viral.report"

#Run Diamond using the CARD database: f1 reads
echo $dmndHeaders > $prefix".f1.card"
diamond blastx --outfmt 6 $dmndHeaders --db $cardDB --evalue $evalue --query $f1 >> $prefix".f1.card" 

#Run Diamond using the CARD database: f2 reads
echo $dmndHeaders > $prefix".f2.card"
diamond blastx --outfmt 6 $dmndHeaders --db $cardDB --evalue $evalue --query $f2 >> $prefix".f2.card"                      

#Run Diamond using the VFDB database: f1 reads
echo $dmndHeaders > $prefix".f1.vfdb"
diamond blastx --outfmt 6 $dmndHeaders --db $vfdbDB --evalue $evalue --query $f1 >> $prefix".f1.vfdb"

#Run Diamond using the VFDB database: f2 reads
echo $dmndHeaders > $prefix".f2.vfdb"
diamond blastx --outfmt 6 $dmndHeaders --db $vfdbDB --evalue $evalue --query $f2 >> $prefix".f2.vfdb"

#run MaxBin
mkdir tmp_dir
$maxbin -thread $threads -contig $outdir/contigs.fasta -reads $f1 -reads2 $f2 -out tmp_dir/maxbin
mv tmp_dir/maxbin.*.fasta $outdir
rm tmp_dir/*
rm -r tmp_dir

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
    if [ ! -s "tmp_mash.txt" ]
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

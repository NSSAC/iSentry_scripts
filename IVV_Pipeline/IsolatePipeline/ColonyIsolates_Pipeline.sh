
#Set path for scripts before running
#Make sure to set the same variable in the SetupMash_Isolates.sh script
rel_path="/scratch/cc8dm/iSentry_scripts/IVV_Pipeline/IsolatePipeline/"

#Assuming reads files do not need to be processed after assembling into contigs

#LOAD MODULES: gcc MUST be loaded before diamond
module load gcc
module load diamond

#Number of threads
threads="1"

#process arguments
sample="none"
card="0"
vfdb="0"
outdir="IsolateOutput"
forward="none"
reverse="none"
while getopts "ts:cvof:r" arg; do
    case $arg in
        t)
            threads=${OPTARG}
            ;;
        s)
            sample=${OPTARG}
            ;;
        c)
            card="1"
            ;;
        v)
            vfdb="1"
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

#check if forward reads are provided, and then if it exists
if [ $forward == "none" ] 
then
    echo "At least one reads file must be provided: [ -f FILE ]" 1>&2
    exit -1
fi
if [ ! -f $forward ]
then
    echo "$forward does not exist" 1>&2
    exit -1
fi

#Create Output directory if it doesn't exist
if [ ! -d $outdir ]
then
    mkdir $outdir
fi
cd $outdir

#Run Spades
sh "$rel_path"RunSpades_Isolates.sh -t $threads -f $forward -r $reverse
#check if contigs file exists
if [ ! -f tmp_dir/contigs.fasta ]
then
    echo "spades failed for isolate $sample" 1>&2
    exit -1
else
    mv tmp_dir/contigs.fasta .
    rm -r tmp_dir
fi 

#Run Checkm script on contigs file
sh "$rel_path"RunCheckm_Isolates.sh

#Run Diamond: card and vfdb
diamond_cmd = "sh ""$rel_path""RunDiamond_Isolates.sh"
if [ $card == "1" ] 
then 
    diamond_cmd=$diamond_cmd" -c" 
fi
if [ $vfdb == "1" ]
then
    diamond_cmd=$diamond_cmd" -v" 
fi
if [ $card == "0" ] && [ $vfdb == "0" ]
then
    echo "Diamond against both card and vfdb has been turned off" 1>&2
else
    diamond_cmd=$diamond_cmd" -i "$sample" -f contigs.fasta" 
    eval $diamond_cmd
    #sh RunDiamond_Isolates.sh -c -v -i $sample -f $forward
fi

#Run Mash
sh "$rel_path"SetupMash_Isolates.sh -t $threads -f contigs.fasta -s $sample 

#Run Kraken script
sh "$rel_path"RunKraken.sh -t $threads -f contigs.fasta -i $sample -b -r -g -s

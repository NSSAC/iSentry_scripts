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
while getopts "t:s:cvo:f:r:" arg; do
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

#TODO: Change relative paths for shells and python scripts to absolute paths 
#TODO: process differently for contigs and reads. Maybe copy and make a different script
#Can assume there will be two reads per file since the others had that as well, change to fit that assumption
#If there are samples with only one reads file then make a separate file for them or use the contigs pipeline

#Create Output directory if it doesn't exist
if [ ! -d $outdir ]
then
    mkdir $outdir
fi
cd $outdir

#Run Spades
sh ../RunSpades_Isolates.sh -t $threads -f $forward -r $reverse
#check if contigs file exists
if [ ! -f tmp_dir/contigs.fasta ]
then
    echo "spades failed for isolate $sample" 1>&2
    exit -1
else
    mv tmp_dir/contigs.fasta .
    rm -r tmp_dir
    cd ..
fi 

#Run Diamond: card and vfdb
#TODO: replace with string building and only have two run commands
if [ $card == "1" ] && [ $vfdb == "1" ]
then 
    sh RunDiamond_Isolates.sh -c -v -i $sample -f $forward
    if [ $reverse != "none" ]
    then
        sh RunDiamond_Isolates.sh -c -v -i $sample -f $reverse
    fi
elif [ $card == "1" ] && [ $vfdb == "0" ]
then

elif [ $card == "0" ] && [ $vfdb == "1" ]
then

else
    echo "Diamond against both card and vfdb has been turned off" 1>&2
fi

#Run Mash
sh SetupMash_Isolates.sh -t $threads -f $foward -s $sample 

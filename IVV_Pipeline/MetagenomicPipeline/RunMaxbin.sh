
maxbin="/project/biocomplexity/isentry/src/MaxBin-2.2.7/run_MaxBin.pl"

#arguments
threads="1"
contigs="none"
forward="none"
reverse="none"
outdir="Bin_Output"
sample="maxbin"
while getopts "t:c:f:r:o:s:" arg; do
    case $arg in
        t)
            threads=$OPTARG
            ;;
        c)
            contigs=$OPTARG
            ;;
        f)
            forward=$OPTARG
            ;;
        r)
            reverse=$OPTARG
            ;;
        o)
            outdir=$OPTARG
            ;;
        s)
            sample=$OPTARG
            ;;
        \?)
            echo "Invalid option: -$OPTARG" 1>&2
            exit -1
            ;;
    esac  
done

#Check for contigs and forward read
if [[ $contigs == "none" ]] || [[ $forward == "none" ]]
then
    echo "Pass in contigs and reads: [ -c FILE ] [ -f FILE ]" 1>&2
    exit -1
fi
if [[ ! -f $contigs ]]
then
    echo "$contigs does not exist" 1>&2
    exit -1
fi
if [[ ! -f $forward ]]
then
    echo "$forward does not exist" 1>&2
    exit -1
fi

mkdir tmp_dir
if [[ "$reverse" == "none" ]] 
then
    $maxbin -thread $threads -contig $contigs -reads $forward  -out tmp_dir/$sample
else
    $maxbin -thread $threads -contig $contigs -reads $forward -reads2 $reverse -out tmp_dir/$sample
fi
mv tmp_dir/$sample.*.fasta $outdir
rm -r tmp_dir

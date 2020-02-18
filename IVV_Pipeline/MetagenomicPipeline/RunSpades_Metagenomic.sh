forward="none"
reverse="none"
threads="1"
while getops "f:rt" arg; do
    case $arg in 
        f)
            forward=${OPTARG}
            ;;
        r)
            reverse=${OPTARG}
            ;;
        t)
            threads=${OPTARG}
            ;;
        \?)
            echo "Invalid option: -$OPTARG" 1>&2
            exit -1
            ;;

    esac
done 

#Assemble Fastq Files using SPAdes, remove all but the assembled contigs file
#Flag: --meta, indicates that the fastq files were generated from fastq reads
spades.py --meta -t $threads -1 $f1 -2 $f2 -o tmp_dir 
mv tmp_dir/contigs.fasta $outdir 
rm -r tmp_dir

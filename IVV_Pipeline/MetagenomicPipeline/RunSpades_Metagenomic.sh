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

#check for reads file and that it exists
if [ $forward == "none" ]
then
    echo "Must provide at least one reads file: [ -f FILE ]" 1>&2
    exit -1
fi
if [ ! -f $forward ]
then
    echo "$forward does not exist" 1>&2
    exit -1
fi

#Flag: --meta, indicates that the fastq files were generated from fastq reads
if [ $reverse == "none" ]
then
    spades.py --meta -t $threads -1 $forward  -o tmp_dir 
elif [ ! -f $reverse ]
then
    echo "$reverse does not exist" 1>&2
    exit -1
else
    spades.py --meta -t $threads -1 $forward -2 $reverse -o tmp_dir 
fi

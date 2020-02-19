
bacteriaDB="/project/biocomplexity/isentry/ref_data/kraken2/bacteria"
rdpDB="/project/biocomplexity/isentry/ref_data/kraken2/rdp/"
ggDB="/project/biocomplexity/isentry/ref_data/kraken2/greengenes/"
silvaDB="/project/biocomplexity/isentry/ref_data/kraken2/silva/"

threads="1"
bacteria="0"
rdp="0"
gg="0"
silva="0"
file="none"
sample="isolate"
while getops "f:brgsti" arg; do
    case $arg in 
        f)
            file=${OPTARG}
            ;;
        b)
            bacteria="1"
            ;;
        r)
            rdp="1"
            ;;
        g)
            gg="1"
            ;;
        s)
            silva="1"
            ;;
        t)
            threads=${OPTARG}
            ;;
        i)
            sample=${OPTARG}
            ;;
        \?)
            echo "RunKraken: Invalid option -$OPTARG" 1>&2 
            exit -1
            ;;
    esac
done

if [ $file == "none" ]
then
    echo "Provide a reads file: [ -f FILE ]" 1>&2
    exit -1
fi
if [ ! -f $file ]
then
    echo "$file does not exist" 1>&2
    exit -1
fi  

if [ $bacteria == "1"] 
then
    kraken2 --db $bacteriaDB --threads $threads $file --output $sample"_Bacteria.output" --report $sample"_Bacteria.report"
fi
if [ $rdp == "1" ]
then
    kraken2 --db $rdpDB --threads $threads $file --output $sample"_RDP.output" --report $sample"_RDP.report"
fi
if [ $gg == "1" ] 
then
    kraken2 --db $ggDB --threads $threads $file --output $sample"_GreenGenes.output" --report $sample"_GreenGenes.report"
fi
if [ $silva == "1" ]
then
    kraken2 --db $silvaDB --threads $threads $file --output $sample"_Silva.output" --report $sample"_Silva.report"
fi

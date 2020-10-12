#Run diamond against a card and vfdb database

#Diamond databases
cardDB="/project/biocomplexity/isentry/ref_data/card/card_protein_variant_DB.dmnd"
vfdbDB="/project/biocomplexity/isentry/ref_data/vfdb/VFDB_protein_DB.dmnd"

#Diamond headers
dmndHeaders="qseqid sseqid  length  evalue  bitscore    stitle"

#Default isolate for files
isolate="isolate"

#Diamond evalue threshold
evalue="0.000001"
file="none"
card="0"
vfdb="0"
while getopts "cvh:e:s:f:" args; do
    case $args in 
        c)
            card="1"
            ;;
        v)
            vfdb="1"
            ;;
        h)
            dmndHeaders=${OPTARG}
            ;;
        e)
            evalue=${OPTARG}
            ;;
        i)
            isolate=${OPTARG}
            ;;
        f)
            file=${OPTARG}
            ;;
        \?)
            echo "Invalid option: -$OPTARG" 1>&2
            exit -1
            ;;
    esac
done
#check if file was specified
if [ $file == "none" ]
then
    echo "Provide a query file: [ -f FILE ]" 1>&2
    exit -1
fi

#check if file exists
if [ ! -f $file ]
then
    echo "$file does not exist: exiting" 1>&2
    exit -1
fi

#check if databases are specified and run if they exist
#TODO: Can be rewritten to support a generic database instead of two specific ones
if [ $cardDB == "1" ]
then
    if [ ! -f $cardDB ]
    then
        echo "$cardDB does not exist. Not running diamond with card database" 1>&2
    else
        echo $dmndHeaders > $isolate".card"
        diamond blastx --outfmt 6 $dmndHeaders --db $cardDB --evalue $evalue --query $file >> $isolate".card"
    fi
fi
if [ $vfdb == "1" ] 
then
    if [ ! -f $vfdbDB ]
    then
        echo "$vfdbDB does not exist. Not running diamond with vfdb database" 1>&2
    else
        echo $dmndHeaders > $isolate".vfdb"
        diamond blastx --outfmt 6 $dmndHeaders --db $vfdbDB --evalue $evalue --query $file >> $isolate".vfdb"
    fi
fi

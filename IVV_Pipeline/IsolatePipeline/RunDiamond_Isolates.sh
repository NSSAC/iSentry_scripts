#Run diamond against a card and vfdb database

#prefix for files
prefix="isolate"

#Card and vfdb database: default none to ignore them. Pass database in calling file to run them
cardDB="none"
vfdb="none"

#Desired headers for diamond output
dmndHeaders="qseqid sseqid  length  evalue  bitscore    stitle"

#Diamond: evalue threshold
evalue="0.000001"
file="none"
while getops "c:v:h:e:i:f:" args: do
    case $arg in 
        c)
            cardDB=${OPTARG}
            ;;
        v)
            vfdbDB=${OPTARG}
            ;;
        h)
            dmndHeaders=${OPTARG}
            ;;
        e)
            evalue=${OPTARG}
            ;;
        i)
            prefix=${OPTARG}
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
if [ $cardDB != "none" ]
then
    if [ ! -f $cardDB ]
    then
        echo "$cardDB does not exist. Not running diamond with card database" 1>&2
    else
        echo $dmndHeaders > $prefix".card"
        diamond blastx --outfmt 6 $dmndHeaders --db $cardDB --evalue $evalue --query $file >> $prefix".card"
    fi
fi
if [ $vfdb != "none" ] 
then
    if [ ! -f $vfdbDB ]
    then
        echo "$vfdbDB does not exist. Not running diamond with vfdb database" 1>&2
    else
        echo $dmndHeaders > $prefix".vfdb"
        diamond blastx --outfmt 6 $dmndHeaders --db $vfdbDB --evalue $evalue --query $file >> $prefix".vfdb"
    fi
fi

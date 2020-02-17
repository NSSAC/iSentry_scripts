#Run mash dist on isolate samples. Can be reads or isolates

#Database paths for Patric (Mash)
patricDB="/project/biocomplexity/isentry/ref_data/mash/patric_all.msh"
patricMapping="patric_genomes_names.txt"

#Headers for diamond and mash output
mashHeaders="gene-ID    species  distance   p-value shared-hashes"

#additional arguments
threads="1"
file="none"
sample="isolate"
while getopts "d:m:h:t:f:s:" arg: do
    case $arg in 
        t)
            threads=${OPTARG}
            ;;
        d)
            patricDB=${OPTARG}
            ;; 
        m)
            patricMapping=${OPTARG}
            ;;
        h)
            mashHeaders=${OPTARG}
            ;;
        f)
            file=${OPTARG}
            ;;
        s)
            sample=${OPTARG}
            ;;
        \?)
            echo "Invalid option: -$OPTARG" 1>&2
            exit -1
    esac
done

#Check that file is specified and exists
if [ $file == "none" ]
then
    echo "Error: pass in a reads file [ -f FILE ]" 1>&2
    exit -1
fi
if [ ! -f $file ]
then
    echo "Error: $file does not exist" 1>&2
    exit -1
fi

#Check that database files exists
if [ ! -f $patricDB ]
then
    echo "patric-mash database does not exist:$patricDB" 1>&2
    exit -1
fi
if [ ! -f $patricMapping ]
then
    echo "patric-mapping file does not exist:$patricMapping" 1>&2
    exit -1
fi

#Run Mash Dist on the colony isolate samples
python RunMash_Isolates.py -t $threads -d $patricDB -f $file -s $sample

#Filter Mash Dist values based on threshold values
python MapIsolatesToPatricGenome.py -f $sample".mash_out.txt" -m $patricMapping -o $sample".PatricMapping.txt" -n 2

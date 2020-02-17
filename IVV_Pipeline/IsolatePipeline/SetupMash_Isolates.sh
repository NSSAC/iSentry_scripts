#Run mash dist on isolate samples. Can be reads or isolates

#Database paths for Patric (Mash)
patricDB="/project/biocomplexity/isentry/ref_data/mash/patric_all.msh"
patricMapping="patric_genomes_names.txt"

#Headers for diamond and mash output
mashHeaders="gene-ID    species  distance   p-value shared-hashes"

#mash threshold: Sets maximum distance value for output a mash pair
mash_threshold="0.05"

#additional arguments
threads="1"
file="none"
sample="isolate"
while getopts "d:m:h:t:f:" arg: do
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
    echo "Error: pass in a reads file [ -f FILE]" 1>&2
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

python RunMash_Isolates.py -t $threads -d $patricDB -m $patricMapping -c $mashHeaders -f $file -s $sample

exit 1
#Run mash: patric
#TODO: Replace this with a function that loops over possible distance values finding the minimum value
#TODO: Create a python script that calls mash as a subprocees and parse through the file using this script
mash dist -d $mash_threshold $patricDB $contigs | cut -f 1,3,4,5  > tmp_mash.txt 
if [ ! -s "tmp_mash.txt" ]
then
    echo "Mash failed at threshold $mash_threshold" > mash_failed.txt
else
    grep -F -w "$(cut -f 1 tmp_mash.txt | cut -d'/' -f 6)" $patricMapping > tmp_patric.txt 
    echo $mashHeaders > $prefix".patric"
    paste tmp_patric.txt tmp_mash.txt | cut -f 1,2,4,5,6 >> $prefix".patric" 
    rm tmp_patric.txt
fi 
rm tmp_mash.txt

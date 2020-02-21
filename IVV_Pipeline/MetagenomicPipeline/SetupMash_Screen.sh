#Run mash dist on isolate samples. Can be reads or isolates

#Also set path
abs_path="/scratch/cc8dm/iSentry_scripts/IVV_Pipeline/MetagenomicPipeline/"

#Database paths for Patric (Mash)
patricDB="/project/biocomplexity/isentry/ref_data/mash/patric_all.msh"
patricMapping="$abs_path""patric_genomes_names.txt"

#Headers for diamond and mash output
mashHeaders="gene-ID    species  distance   p-value shared-hashes"

#additional arguments
threads="1"
file="none"
sample="isolate"
while getopts "d:m:t:f:s:" arg; do
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
if [[ $file == "none" ]]
then
    echo "Error: pass in a reads file [[ -f FILE ]]" 1>&2
    exit -1
fi
if [[ ! -f $file ]]
then
    echo "Error: $file does not exist" 1>&2
    exit -1
fi

#Check that database files exists
if [[ ! -f $patricDB ]]
then
    echo "patric-mash database does not exist:$patricDB" 1>&2
    exit -1
fi
if [[ ! -f $patricMapping ]]
then
    echo "patric-mapping file does not exist:$patricMapping" 1>&2
    exit -1
fi

#Run Mash Screen on metagenomic sample
python "$abs_path"RunMash_Screen.py -t $threads -d $patricDB -f $file -s $sample

#Check that return value indicates the program completed successfully: exit == 0
ret_value=$?
if [[ $ret -ne 0 ]]
then
    echo "Error occurred during Mash screen for $file" 1>&2
    exit 1
fi

#Filter Mash Screen values based on threshold values
python "$abs_path"MapIsolatesToPatricGenome.py -f $sample".mash_screen.out" -m $patricMapping -o $sample".PatricMapping.txt"

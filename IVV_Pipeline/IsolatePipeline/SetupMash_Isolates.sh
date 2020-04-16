#Run mash dist on isolate samples. Can be reads or isolates

#Also set path
rel_path="/scratch/cc8dm/iSentry_scripts/IVV_Pipeline/IsolatePipeline/"

#Database paths for Patric (Mash)
patricDB="/project/biocomplexity/isentry/ref_data/mash/patric_all.msh"
patricMapping="$rel_path""patric_genomes_names.txt"

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

#Run Mash Dist on the colony isolate samples
python "$rel_path"RunMash_Isolates.py -t $threads -d $patricDB -f $file -s $sample

#Filter Mash Dist values based on threshold values
python "$rel_path"MapIsolatesToPatricGenome.py -f $sample".mash_dist.out" -m $patricMapping -o $sample".PatricMapping.txt" -n 2

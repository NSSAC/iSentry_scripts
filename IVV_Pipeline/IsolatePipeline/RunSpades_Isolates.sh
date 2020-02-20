#Runs spades on colony isolates and moves pertinent files into the specified directory

#Set arguments
threads=1
file1="none"
file2="none"
while getopts "t:f:r:" arg; do
    case $arg in 
        t)  
            threads=${OPTARG}
            ;; 
        f)
            file1=${OPTARG}
            ;;
        r) 
            file2=${OPTARG}
            ;;
        \?)
            echo "Invalid option: -$OPTARG" 1>&2
            exit -1
            ;;
    esac
done

#Check for reads file and that it exists
if [[ $file1 == "none" ]] 
then
    echo "Must provide at least one reads file: [ -f FILE ]" 1>&2
    exit -1
fi
if [[ ! -f $file1 ]]
then
    echo "$file1 does not exist" 1>&2
    exit -1
fi

#Run single end spades if one file. Else, run paired end spades
if [[ $file2 == "none" ]]
then
    spades.py -t $threads -f $file1 -o tmp_dir
elif [[ ! -f $file2 ]]
then
    echo "$file2 does not exist" 1>&2
    exit -1
else
    spades.py -t $threads -1 $file1 -2 $file2 -o tmp_dir
fi

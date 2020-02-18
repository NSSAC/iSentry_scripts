#Runs spades on colony isolates and moves pertinent files into the specified directory

#Set arguments
threads=1
file1="none"
file2="none"
while getopts "tf:r" arg; do
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

#Check for reads file
if [ $file1 == "none" ] 
then
    echo "Must provide at least one reads file: [ -f FILE ]" 1>&2
    exit -1
fi

#Run single end spades if one file. Else, run paired end spades
if [ $file2 == "none" ]
then
    spades.py -t $threads -f $file1 -o tmp_dir
else
    spades.py -t $threads -1 $file1 -2 $file2 -o tmp_dir
fi

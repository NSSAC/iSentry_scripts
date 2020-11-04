#!/bin/bash

#SBATCH -n 1
#SBATCH --time=168:00:00
#SBATCH -p bii
#SBATCH -A isentry
#SBATCH --mem=150000

#SBATCH -o Bracken.out
#SBATCH -e Brakcen.err

# modules
module load gcc

bacteriaDB="/project/biocomplexity/isentry/ref_data/kraken2/bacteria/"

/project/biocomplexity/isentry/src/Bracken/bracken -d $bacteriaDB -i $1 -o $2 -r 251 -l S 

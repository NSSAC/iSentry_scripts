import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f','--mashfile')
parser.add_argument('-m','--mapfile')
parser.add_argument('-o','--outfile')
#parser.add_argument('-n','--num_thresholds',default=1)

args = parser.parse_args()

#Read in genomes from mapping file
patric_genomes = args.mapfile
patric_dict = {}
with open(patric_genomes,"r") as pg:
    next(pg)
    for line in pg:
        line = line.rstrip().split("\t")
        patric_dict[line[0]] = line[1]

#Read in mash file and map to patric mapping file
infile = args.mashfile
data_list = []
with open(infile,"r") as f:
    for line in f:
        line = line.rstrip().split()
        genome_id = line[4].split("/")[-2]
        if genome_id in patric_dict:
            data = (patric_dict[genome_id],genome_id,line[0],line[1],line[2],line[3])
            data_list.append(data)

#Output info based on sorted values
outfile = args.outfile
with open(outfile,"w") as o:
    o.write("Genome\tGenomd_ID\tContainment\tShared-Hashes\tMedian-Multiplicity\tP-Value\n")
    for data in data_list:
        o.write("%s\n"%("\t".join(data)))

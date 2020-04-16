import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-p','--path')
parser.add_argument('-o','--outfile')

args = parser.parse_args()

data_list = []

for genus_file in glob.glob(args.path):
    isolate = genus_file.split("/")[0]
    with open(genus_file,"r") as gf:
        for line in gf:
            line = line.rstrip().split()
            if len(line) == 7:
                line[5] = line[5]+"-"+line[6]
                line.pop()
            if float(line[0]) > 0.0:
                line.append(isolate)
                data_list.append(line)

with open(args.outfile,"w") as o:
    o.write("Isolate\tMapped_Pct\tReads_Covered\tReads_Assigned\tRank\tTaxon_ID\tGenus\n")
    for data in data_list:
        o.write("%s\t%s\n"%(data[-1],"\t".join(data[:len(data)-1])))

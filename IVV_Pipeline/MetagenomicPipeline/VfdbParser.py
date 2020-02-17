import argparse
import sys

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-g","--id",help="Gene id")
parser.add_argument("-f","--fof",help="Friend or Foe")
parser.add_argument("-v","--vfdb",help="VFDB input File")
parser.add_argument("-o","--output",help="Output File Name",default="vfdboutput.txt")

args = parser.parse_args()

geneid = args.id
friendOrFoe = args.fof
vfdbfile = args.vfdb
outfile = args.output

qseqid,sseqid,length,evalue,bitscore,stitle=0,1,2,3,4,5

gene_dict = {}
with open(vfdbfile,"r") as v:
    for line in v:
        gene_info = line.rstrip().split() 
        gene = gene_info[6][1:-1]
        if gene not in gene_dict:
            gene_dict[gene] = 1
        else:
            gene_dict[gene] = gene_dict[gene] + 1

with open(outfile,"a") as o:
    for gene in gene_dict:
        tmp = gene+"\t"+str(gene_dict[gene])
        o.write("%s\t%s\tvfdb\t%s\n"%(geneid,friendOrFoe,tmp))

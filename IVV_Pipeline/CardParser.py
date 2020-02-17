import argparse
import sys

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-g","--id",help="Gene id")
parser.add_argument("-f","--fof",help="Friend or Foe")
parser.add_argument("-c","--card",help="Card Input File")
parser.add_argument("-o","--output",help="Output File Name",default="cardoutput.txt")

args = parser.parse_args()

geneid = args.id
friendOrFoe = args.fof
cardfile = args.card
outfile = args.output

qseqid,sseqid,length,evalue,bitscore,stitle=0,1,2,3,4,5

gene_dict = {}
with open(cardfile,"r") as c:
    for line in c:
        line = " ".join(line.rstrip().split()).split(" ",5)
        title_info = line[stitle].split("|")[-1]
        #1: test for line with just the gene name
        if title_info.split(" ",1)[1][0] == "[":
            gene = title_info.split(" ",1)[0]
            if gene not in gene_dict:
                gene_dict[gene] = {}
                gene_dict[gene]["drug"] = "NA"
                gene_dict[gene]["count"] = 1
            else:
                gene_dict[gene]["count"] = gene_dict[gene]["count"] + 1
        elif title_info.split(" ",1)[0] == "antibiotic":
            gene = title_info.split(" ")[2]
            if gene not in gene_dict:
                gene_dict[gene] = {}
                gene_dict[gene]["drug"] = "NA"
                gene_dict[gene]["count"] = 1
            else:
                gene_dict[gene]["count"] = gene_dict[gene]["count"] + 1
        else:
            info = title_info.split("[")[0].rstrip().split(" ")
            if len(info) <= 2:
                continue
            gene = info[2]
            if gene == "mutant":
                continue
            if gene == "intrinsic":
                if gene not in gene_dict:
                    gene_dict[info[3]] = {}
                    gene_dict[info[3]]["drug"] = info[-1]
                    gene_dict[info[3]]["count"] = 1
                else:
                    gene_dict[info[3]]["count"] = gene_dict[info[3]]["count"] + 1
            if info[-1] == "resistance":
                if info[-2] == "antibiotic":
                    if info[-3] == "conferring":
                        continue
                    if gene not in gene_dict:
                        gene_dict[gene] = {}
                        gene_dict[gene]["drug"] = info[-3]
                        gene_dict[gene]["count"] = 1
                    else:
                        gene_dict[gene]["count"] = gene_dict[gene]["count"] + 1
                else:
                    if gene not in gene_dict:
                        gene_dict[gene] = {}
                        gene_dict[gene]["drug"] = info[-2]
                        gene_dict[gene]["count"] = 1
                    else:
                        gene_dict[gene]["count"] = gene_dict[gene]["count"] + 1
            else:
                if info[-1] == "acid" or info[-1] == "antibiotics":
                    if gene not in gene_dict:
                        gene_dict[gene] = {}
                        gene_dict[gene]["drug"] = " ".join(info[-2:])
                        gene_dict[gene]["count"] = 1
                    else:
                        gene_dict[gene]["count"] = gene_dict[gene]["count"] + 1
                else:
                    if info[-1] == "mutation" or info[-1] == "mutants":
                        continue
                    if gene not in gene_dict:
                        gene_dict[gene] = {}
                        gene_dict[gene]["drug"] = info[-1]
                        gene_dict[gene]["count"] = 1
                    else:
                        gene_dict[gene]["count"] = gene_dict[gene]["count"] + 1

with open(outfile,"a") as o:
    for gene in gene_dict:
        tmp = gene+"\t"+gene_dict[gene]["drug"]+"\t"+str(gene_dict[gene]["count"])
        o.write("%s\t%s\tcard\t%s\n"%(geneid,friendOrFoe,tmp))

import sys
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-f","--file",help="Counts file created from CardParser.py and VfdbParser.py")
parser.add_argument("-c","--card",help="Card Output Matrix filename",default="card_matrix.txt")
parser.add_argument("-v","--vfdb",help="VFDB Output Matrix filename",default="vfdb_matrix.txt")

args = parser.parse_args()

infoFile = args.file
vfdbMatrix = args.vfdb
cardMatrix = args.card

cardDict = {}
vfdbDict = {}
all_vfdb_genes = []
all_card_pairs = []
all_org_ids = []
with open(infoFile,"r") as info:
    for line in info:
       line = line.rstrip().split("\t")
       if line[0] not in all_org_ids:
           all_org_ids.append(line[0])
       if line[2] == "vfdb":
           if line[3] not in all_vfdb_genes:
               all_vfdb_genes.append(line[3])
           if line[0] not in vfdbDict:
               vfdbDict[line[0]] = {}
           vfdbDict[line[0]][line[3]] = line[4]
       else: #card
           g_d = line[3]+"_"+line[4]
           if g_d not in all_card_pairs:
               all_card_pairs.append(g_d)
           if line[0] not in cardDict:
               cardDict[line[0]] = {}
           cardDict[line[0]][g_d] = line[5]

vfdbZeros = ["0"]*len(all_vfdb_genes)
cardZeros = ["0"]*len(all_card_pairs)

with open(vfdbMatrix,"w") as vM:
    #write headers
    vM.write("Organism ID")
    for gene in all_vfdb_genes:
        vM.write("\t%s"%gene)
    vM.write("\n")
    #write counts
    for ID in all_org_ids:
        vM.write("%s"%ID)
        if ID not in vfdbDict:
            vM.write("\t%s"%"\t".join(vfdbZeros))
        else:
            for gene in all_vfdb_genes:
                if gene not in vfdbDict[ID]:
                    vM.write("\t0")
                else:
                    vM.write("\t%s"%vfdbDict[ID][gene])
        vM.write("\n")

with open(cardMatrix,"w") as cM:
    #write headers
    cM.write("Organism ID")
    for gene in all_card_pairs:
        cM.write("\t%s"%gene)
    cM.write("\n")
    #write counts
    for ID in all_org_ids:
        cM.write("%s"%ID)
        if ID not in cardDict:
            cM.write("\t%s"%"\t".join(cardZeros))
        else:
            for gene in all_card_pairs:
                if gene not in cardDict[ID]:
                    cM.write("\t0")
                else:
                    cM.write("\t%s"%cardDict[ID][gene])
        cM.write("\n")

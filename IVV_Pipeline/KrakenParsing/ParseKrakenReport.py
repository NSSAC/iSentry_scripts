import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-r','--report')

args = parser.parse_args()
rank_index = 3
genus_list = [] 
with open(args.report,"r") as r:
    for line in r:
        line = line.rstrip().split()
        info = [s.rstrip() for s in line]
        if info[rank_index] == "G":
            if len(info) > 6:
                print(line)
                info[5] == info[5] + " " + info[6]
                list.pop(info)
            genus_list.append(info)

#print(genus_list)

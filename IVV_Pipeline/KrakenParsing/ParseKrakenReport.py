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
                info[5] = info[5] + " " + info[6] 
                list.pop(info)
            info[0] = float(info[0])
            info = tuple(info)
            genus_list.append(info)

genus_list_sorted = sorted(genus_list,key = lambda x: x[0],reverse=True)



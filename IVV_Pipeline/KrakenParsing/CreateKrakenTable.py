import glob

path = "*/*GenusOutput.txt"

data_list = []

for genus_file in glob.glob(path):
    isolate = genus_file.split("/")[0]
    isolate = isolate.split("_")[1]
    with open(genus_file,"r") as gf:
        for line in gf:
            line = line.rstrip().split()
            if len(line) == 7:
                line[5] = line[5]+"-"+line[6]
                line.pop()
            if float(line[0]) > 0.0:
                line.append(isolate)
                data_list.append(line)

outfile = "IVV_August2019_16s_Genus_Table.txt"
with open(outfile,"w") as o:
    o.write("Isolate\tMapped_Pct\tReads_Covered\tReads_Assigned\tRank\tTaxon_ID\tGenus\n")
    for data in data_list:
        o.write("%s\t%s\n"%(data[-1],"\t".join(data[:len(data)-1])))

#!/home/cc8dm/miniconda3/bin/python

import os,shutil,subprocess,sys,argparse,time

parser = argparse.ArgumentParser()
parser.add_argument("-r1","--read1",required=True)
parser.add_argument("-r2","--read2",required=True)
parser.add_argument("-t","--threads",default="1")
parser.add_argument("-o","--outdir",default="./")
parser.add_argument("-p","--prefix",required=True)

args = parser.parse_args()
print(args)

#Mash database and id mapping
patricDB = "/project/biocomplexity/isentry/ref_data/mash/patric_all.msh"
patricMapping="/project/biocomplexity/isentry/ref_data/mash/patric_genomes_names.txt"

#Headers for mash
mashHeaders=["gene-ID","species","distance","p-value"]
mash_threshold = 0.25

#Check programs are loaded
if not shutil.which("diamond"):
    sys.exit("diamond not in path")
if not shutil.which("kraken2"):
    sys.exit("kraken not in path")

#Setup time collection dictionary
time_programs = ["Mash","Kraken","Diamond","Spades","Maxbin"]
time_dict = {}
for program in time_programs:
    time_dict[program] = {}
    if program != "Spades":
        time_dict[program]["Contigs"] = {}
        time_dict[program]["Bins"] = {}

#TODO: I don't think these should be run on the reads but keep for now
#RunKraken on forward reads
#Run Mash on forward reads
#Run Diamond on forward reads

#RunSpades_Metagenomic.sh
spades_dir = args.prefix+"_spades"
if not os.path.isdir(spades_dir):
    os.mkdir(spades_dir)
#os.mkdir(spades_dir)
contig_file = os.path.join(spades_dir,args.prefix+".fa")
if not os.path.exists(contig_file):
    time_dict["Spades"]["Start"] = time.time()
    subprocess.check_call(["RunSpades_Metagenomic.py","-t",args.threads,"-r1",args.read1,"-r2",args.read2,"-s",spades_dir,"-p",args.prefix])
    time_dict["Spades"]["End"] = time.time()
else:
    print("contigs file exists: skipping spades assembly")
    time_dict["Spades"]["Start"] = "NA"
    time_dict["Spades"]["End"] = "NA"
#Check for existence of contigs file
if not os.path.exists(contig_file):
    sys.stderr.write("Contigs file doesn't exist for sample %s\n"%args.prefix)
    sys.exit()

#Run Kraken on contigs
if False:
    kraken_log_stdout = "kraken.stdout"
    kraken_log_stderr = "kraken.stderr"
    time_dict["Kraken"]["Contigs"]["Start"] = time.time()
    with open(kraken_log_stdout,"w") as klo, open(kraken_log_stderr,"w") as kle:
        subprocess.check_call(["RunKraken.py","-t",args.threads,"-p",args.prefix,"-r",contig_file],stdout=klo,stderr=kle)
    time_dict["Kraken"]["Contigs"]["End"] = time.time()
else:
    print("Skipping Kraken")

#Run Mash on contigs

#Run Diamond on contigs

#Run maxbin
#Within maxbin: run diamond, kraken, and mash


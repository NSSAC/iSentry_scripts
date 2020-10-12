import os,shutil,subprocess,sys,argparse,time

parser = argparse.ArgumentParser()
parser.add_argument("-r1","--read1",required=True)
parser.add_argument("-r2","--read2",required=True)
parser.add_argument("-t","--threads",default="1")
parser.add_argument("-o","--outdir",default="./")
parser.add_argument("-p","--prefix",required=True)

args = parser.parse_args()

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
os.mkdir(spades_dir)
time_dict["Spades"]["Start"] = time.time()
subprocess.check_call(["python","RunSpades_Metagenomic.py","-t",args.threads,"-r1",args.r1,"-r2",args.r2,"-s",spades_dir,"-p",args.prefix])
time_dict["Spades"]["End"] = time.time()
#Check for existence of contigs file
contig_file = os.path.join(spades_dir,"contigs.fasta")
if not os.path.exits(contig_file):
    sys.stderr.write("Contigs file doesn't exist for sample %s\n"%args.prefix)

#Run Kraken on contigs
#Run Mash on contigs
#Run Diamond on contigs

#Run maxbin
#Within maxbin: run diamond, kraken, and mash


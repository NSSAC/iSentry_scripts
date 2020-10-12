import argparse
import subprocess
import time
import os

parser = argparse.ArgumentParser()

parser.add_argument("-t","--threads",default=1,type=str)
parser.add_argument("-p","--prefix",default="sample")
parser.add_argument("-o","--outdir",default="./")
parser.add_argument("-f","--forward",required=True)
parser.add_argument("-r","--reverse",required=True)

args = parser.parse_args()

#Databases
bacteriaDB = "/project/biocomplexity/isentry/ref_data/kraken2/bacteria"
vfdbDB="/project/biocomplexity/isentry/ref_data/vfdb/VFDB_protein_DB.dmnd"
patricDB="/project/biocomplexity/isentry/ref_data/mash/patric_all.msh"
patricMapping = "/scratch/cc8dm/iSentry_scripts/IVV_Pipeline/IsolatePipeline/old_patric_genomes_names.txt"

#Set up time collection data object
time_dict = {}
time_dict["Kraken"] = {}
time_dict["Kraken"]["Reads"] = {}
time_dict["Kraken"]["Contigs"] = {}
time_dict["Diamond"] = {}
time_dict["Diamond"]["Reads"] = {}
time_dict["Diamond"]["Contigs"] = {}
time_dict["PathoProfile"] = {}
time_dict["PathoProfile"]["Reads"] = {}
time_dict["PathoProfile"]["Contigs"] = {}
time_dict["Spades"] = {}
time_dict["Mash"] = {}
time_dict["Mash"]["Reads"] = {}
time_dict["Mash"]["Contigs"] = {}

#Run kraken on reads (forward)
kraken_output_reads = args.prefix+"_reads_Bacteria.output"
kraken_report_reads = args.prefix+"_reads_Bacteria.report"
time_dict["Kraken"]["Reads"]["Start"] = time.time()
with open("Kraken_Reads_Log.txt","w") as log:
    subprocess.check_call(["kraken2","--db",bacteriaDB,"--threads",args.threads,args.forward,"--output",kraken_output_reads,"--report",kraken_report_reads],stdout=log)
time_dict["Kraken"]["Reads"]["End"] = time.time()

#Run Mash on reads (forward)
mash_reads = args.prefix+"_reads"
map_reads = args.prefix+"_reads.PatricMapping.txt"
time_dict["Mash"]["Reads"]["Start"] = time.time()
subprocess.check_call(["python","/scratch/cc8dm/iSentry_scripts/IVV_Pipeline/IsolatePipeline/RunMash_Isolates.py","-t",args.threads,"-d",patricDB,"-f",args.forward,"-s",mash_reads])
time_dict["Mash"]["Reads"]["End"] = time.time()
subprocess.check_call(["python","/scratch/cc8dm/iSentry_scripts/IVV_Pipeline/IsolatePipeline/MapIsolatesToPatricGenome.py","-f",mash_reads+".mash_dist.out","-m",patricMapping,"-o",map_reads,"-n","2"])

#Run diamond on forward reads
evalue = "0.000001"
dmndHeaders = ["qseqid","sseqid","length","evalue","bitscore" ,"stitle"]
vfdb_reads_output = args.prefix+"_contigs.vfdb"
time_dict["Diamond"]["Reads"]["Start"] = time.time()
with open(vfdb_reads_output,"w") as o:
    subprocess.check_call(["diamond","blastx","--outfmt","6"]+dmndHeaders+["--db",vfdbDB,"--evalue",evalue,"--query",args.forward],stdout=o)
time_dict["Diamond"]["Reads"]["End"] = time.time()

#Run pathogenicity profile on reads (forward)
pp_reads_output = args.prefix+"_reads.pp"
time_dict["PathoProfile"]["Reads"]["Start"] = time.time()
with open(pp_reads_output,"w") as o:
    subprocess.check_call(["python","/scratch/cc8dm/iSentry_scripts/IVV_Pipeline/PathogenicityProfile/scoreGenomeByVfdbLogodds.py",vfdb_reads_output],stdout=o)
time_dict["PathoProfile"]["Reads"]["End"] = time.time()

#Run SPAdes
time_dict["Spades"]["Start"] = time.time()
os.mkdir("Contigs")
subprocess.check_call(["python","/scratch/cc8dm/iSentry_scripts/IVV_Pipeline/SPAdes/RunSpades.py","-t",args.threads,"-r",args.forward+","+args.reverse,"-o","Contigs","-p",args.prefix])
contigs_file = os.path.join("Contigs",args.prefix+".fa")
if not os.path.exists(contigs_file):
    sys.stderr.write("Contigs file does not exist")
    sys.exit(-1)
time_dict["Spades"]["End"] = time.time()

#Run Kraken on contigs
kraken_output_contigs = args.prefix+"_contigs_Bacteria.output"
kraken_report_contigs = args.prefix+"_contigs_Bacteria.report"
time_dict["Kraken"]["Contigs"]["Start"] = time.time()
with open("Kraken_Contigs_Log.txt","w") as log:
    subprocess.check_call(["kraken2","--db",bacteriaDB,"--threads",args.threads,contigs_file,"--output",kraken_output_contigs,"--report",kraken_report_contigs],stdout=log)
time_dict["Kraken"]["Contigs"]["End"] = time.time()

#Run Mash on contigs
mash_contigs = args.prefix+"_contigs"
map_contigs = args.prefix+"_contigs.PatricMapping.txt"
time_dict["Mash"]["Contigs"]["Start"] = time.time()
subprocess.check_call(["python","/scratch/cc8dm/iSentry_scripts/IVV_Pipeline/IsolatePipeline/RunMash_Isolates.py","-t",args.threads,"-d",patricDB,"-f",contigs_file,"-s",mash_contigs])
time_dict["Mash"]["Contigs"]["End"] = time.time()
subprocess.check_call(["python","/scratch/cc8dm/iSentry_scripts/IVV_Pipeline/IsolatePipeline/MapIsolatesToPatricGenome.py","-f",mash_contigs+".mash_dist.out","-m",patricMapping,"-o",map_contigs,"-n","2"])

#Run diamond on forward reads
evalue = "0.000001"
dmndHeaders = ["qseqid","sseqid","length","evalue","bitscore" ,"stitle"]
vfdb_contigs_output = args.prefix+"_contigs.vfdb"
time_dict["Diamond"]["Contigs"]["Start"] = time.time()
with open(vfdb_reads_output,"w") as o:
    subprocess.check_call(["diamond","blastx","--outfmt","6"]+dmndHeaders+["--db",vfdbDB,"--evalue",evalue,"--query",contigs_file],stdout=o)
time_dict["Diamond"]["Contigs"]["End"] = time.time()

#Run pathogenicity profile on contigs
pp_contigs_output = args.prefix+"_contigs.pp"
time_dict["PathoProfile"]["Contigs"]["Start"] = time.time()
with open(pp_contigs_output,"w") as o:
    subprocess.check_call(["python","/scratch/cc8dm/iSentry_scripts/IVV_Pipeline/PathogenicityProfile/scoreGenomeByVfdbLogodds.py",vfdb_contigs_output],stdout=o)
time_dict["PathoProfile"]["Contigs"]["End"] = time.time()

#Write out time data
time_file = args.prefix+"_Time.txt"
with open(time_file,"w") as o:
    o.write("Program\tReads\tContigs\n")
    o.write("Kraken\t%s\t%s\n"%(str(time_dict["Kraken"]["Reads"]["End"]-time_dict["Kraken"]["Reads"]["Start"]),str(time_dict["Kraken"]["Contigs"]["End"]-time_dict["Kraken"]["Contigs"]["Start"])))
    o.write("Mash\t%s\t%s\n"%(str(time_dict["Mash"]["Reads"]["End"]-time_dict["Mash"]["Reads"]["Start"]),str(time_dict["Mash"]["Contigs"]["End"]-time_dict["Mash"]["Contigs"]["Start"])))
    o.write("Diamond\t%s\t%s\n"%(str(time_dict["Diamond"]["Reads"]["End"]-time_dict["Diamond"]["Reads"]["Start"]),str(time_dict["Diamond"]["Contigs"]["End"]-time_dict["Diamond"]["Contigs"]["Start"])))
    o.write("PathoProfile\t%s\t%s\n"%(str(time_dict["PathoProfile"]["Reads"]["End"]-time_dict["PathoProfile"]["Reads"]["Start"]),str(time_dict["PathoProfile"]["Contigs"]["End"]-time_dict["PathoProfile"]["Contigs"]["Start"])))
    o.write("Spades\t%s\t-\n"%(str(time_dict["Spades"]["End"]-time_dict["Spades"]["Start"])))

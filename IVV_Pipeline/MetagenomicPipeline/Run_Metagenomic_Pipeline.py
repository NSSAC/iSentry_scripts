#!/home/cc8dm/miniconda3/bin/python

import os,json,shutil,subprocess,sys,argparse,time,glob

parser = argparse.ArgumentParser()
parser.add_argument("-r1","--read1",required=True)
parser.add_argument("-r2","--read2",required=True)
parser.add_argument("-t","--threads",default="1")
parser.add_argument("-o","--outdir",default="./")
parser.add_argument("-p","--prefix",required=True)

args = parser.parse_args()
#print(args)

#Mash database and id mapping
patricDB = "/project/biocomplexity/isentry/ref_data/mash/patric_all.msh"
patricMapping="/project/biocomplexity/isentry/ref_data/mash/patric_genome_name.txt"

#Headers for mash
mashHeaders=["gene-ID","species","distance","p-value"]
mash_threshold = 0.25

#Diamond database and headers
vfdbDB="/project/biocomplexity/isentry/ref_data/vfdb/VFDB_setB_pro.dmnd"
dmndHeaders="qseqid sseqid  length  evalue  bitscore    stitle"

#Check programs are loaded
if not shutil.which("diamond"):
    sys.exit("diamond not in path")
if not shutil.which("kraken2"):
    sys.exit("kraken not in path")

#Setup time collection dictionary
time_programs = ["Mash","Kraken","Bracken","Diamond","Spades","Maxbin","Checkm","PP_Profile"]
time_dict = {}
for program in time_programs:
    time_dict[program] = {}
    if program != "Spades" and program != "Maxbin":
        time_dict[program]["Contigs"] = {}
        time_dict[program]["Bin"] = {}
    if program == "Kraken" or program == "Bracken":
        time_dict[program]["Reads"] = {}

#TODO: RunKraken and bracken on forward reads
kraken_read_log_stdout = args.prefix+".reads.kraken.stdout"
kraken_read_log_stderr = args.prefix+".reads.kraken.stderr"
time_dict["Kraken"]["Reads"]["Start"] = time.time()
kraken_reads_prefix = args.prefix+".reads"
if not os.path.exists(kraken_read_log_stdout):
    print("Running kraken: reads")
    with open(kraken_read_log_stdout,"w") as klo, open(kraken_read_log_stderr,"w") as kle:
        subprocess.check_call(["RunKraken.py","-t",args.threads,"-p",kraken_reads_prefix,"-r",args.read1],stdout=klo,stderr=kle)
    time_dict["Kraken"]["Reads"]["End"] = time.time()

#Run bracken on forward reads
#bracken_read_log_stdout = ".reads.bracken.stdout"
#bracken_read_log_stderr = ".reads.bracken.stderr"
#bracken_bacteriaDB="/project/biocomplexity/isentry/ref_data/kraken2/bacteria/"
#bracken_reads_prefix = args.prefix+".reads"
#time_dict["Kraken"]["Reads"]["Start"] = time.time()
#bracken -d MY_DB -i INPUT -o OUTPUT -r READ_LEN -l LEVEL -t THRESHOLD
#with open(bracken_read_log_stdout,"w") as bro, open(bracken_read_log_stderr,"w") as bre:
#    subprocess.check_call(["bracken","-d",bracken_bacteriaDB,"-i",

#RunSpades_Metagenomic.sh
spades_dir = args.prefix+"_spades"
if not os.path.isdir(spades_dir):
    os.mkdir(spades_dir)
contig_file = os.path.join(spades_dir,args.prefix+".fa")
if not os.path.exists(contig_file):
    print("running spades")
    time_dict["Spades"]["Start"] = time.time()
    subprocess.check_call(["RunSpades_Metagenomic.py","-t",args.threads,"-r1",args.read1,"-r2",args.read2,"-s",spades_dir,"-p",args.prefix])
    time_dict["Spades"]["End"] = time.time()
else:
    print("contigs file exists: skipping spades assembly")
#Check for existence of contigs file
if not os.path.exists(contig_file):
    sys.stderr.write("Contigs file doesn't exist for sample %s\n"%args.prefix)
    sys.exit()

#Run Kraken and Bracken on contigs
if True:
    kraken_log_stdout = "kraken.stdout"
    kraken_log_stderr = "kraken.stderr"
    if not os.path.exists(kraken_log_stdout):
        print("Running kraken: contigs")
        time_dict["Kraken"]["Contigs"]["Start"] = time.time()
        with open(kraken_log_stdout,"w") as klo, open(kraken_log_stderr,"w") as kle:
            subprocess.check_call(["RunKraken.py","-t",args.threads,"-p",args.prefix,"-r",contig_file],stdout=klo,stderr=kle)
        time_dict["Kraken"]["Contigs"]["End"] = time.time()
else:
    print("Skipping Kraken")

#Run Mash Screen on contigs
if True:
    mash_screen_cmd = ["Run_MashScreen.py","-t",args.threads,"-s",args.prefix,"-f",contig_file,"-d",patricDB]
    mash_file = args.prefix+".mash_screen.out"
    if not os.path.exists(mash_file):
        print("running mash: contigs")
        time_dict["Mash"]["Contigs"]["Start"] = time.time()
        subprocess.check_call(mash_screen_cmd)
        time_dict["Mash"]["Contigs"]["End"] = time.time()
else:
    print("Skipping mash")

#Mash mapping
if True:
    mash_mapping_contigs = args.prefix+".mash_screen.mapping" 
    if not os.path.exists(mash_mapping_contigs):
        print("mapping mash contigs")
        mash_screen_map_cmd = ["MapPatricGenomeIDs.py","-m",patricMapping,"-f",mash_file,"-o",mash_mapping_contigs]
        subprocess.check_call(mash_screen_map_cmd)

#Run Diamond on contigs
if True:
    dmnd_contig_output = args.prefix+".diamond.out"
    if not os.path.exists(dmnd_contig_output):
        dmnd_contigs_cmd = ["Run_Diamond.py","-t",args.threads,"-f",contig_file,"-d",vfdbDB,"-p",args.prefix]
        print(dmnd_contigs_cmd)
        time_dict["Diamond"]["Contigs"]["Start"] = time.time()
        subprocess.check_call(dmnd_contigs_cmd)
        time_dict["Diamond"]["Contigs"]["End"] = time.time()

#Checkm
#checkm_cmd = ["/project/biocomplexity/isentry/labkey_env/bin/checkm","-f","./checkm_results","-x","fasta","--individual_markers","domain","Bacteria","./","./Checkm_Results"]
#time_dict["Checkm"]["Contigs"]["Start"] = time.time()
#subprocess.check_call(checkm_cmd)
#time_dict["Checkm"]["Contigs"]["End"] = time.time()

#Run maxbin
if True:
    if not os.path.exists(args.prefix+".001.fasta"):
        maxbin_cmd = ["Run_Maxbin.py","-t",args.threads,"-r1",args.read1,"-r2",args.read2,"-c",contig_file,"-p",args.prefix]
        print(maxbin_cmd)
        time_dict["Maxbin"]["Start"] = time.time()
        subprocess.check_call(maxbin_cmd)
        time_dict["Maxbin"]["End"] = time.time()

#For each bin: run diamond, kraken, and mash
for bin_file in glob.glob(args.prefix+".*.fasta"):
    print("Processing %s"%bin_file)
    bin_prefix = bin_file.replace(".fasta","")
    if not os.path.isdir(bin_prefix):
        os.mkdir(bin_prefix)
    os.chdir(bin_prefix)
    #add to time dictionary
    time_dict["Kraken"]["Bin"][bin_prefix] = {}
    time_dict["Mash"]["Bin"][bin_prefix] = {}
    time_dict["Diamond"]["Bin"][bin_prefix] = {}
    time_dict["PP_Profile"]["Bin"][bin_prefix] = {}
    #kraken and bracken
    kraken_log_stdout = bin_prefix+".kraken.stdout"
    kraken_log_stderr = bin_prefix+".kraken.stderr"
    time_dict["Kraken"]["Bin"][bin_prefix]["Start"] = time.time()
    with open(kraken_log_stdout,"w") as klo, open(kraken_log_stderr,"w") as kle:
        subprocess.check_call(["RunKraken.py","-t",args.threads,"-p",args.prefix,"-r","../"+bin_file],stdout=klo,stderr=kle)
    time_dict["Kraken"]["Bin"][bin_prefix]["End"] = time.time()
    #mash dist
    mash_dist_cmd = ["Run_MashDist.py","-t",args.threads,"-s",bin_prefix,"-f","../"+bin_file,"-d",patricDB]
    time_dict["Mash"]["Bin"][bin_prefix]["Start"] = time.time()
    subprocess.check_call(mash_dist_cmd)
    time_dict["Mash"]["Bin"][bin_prefix]["End"] = time.time()
    bin_mash_file = bin_prefix+".mash_dist.out"
    #mash mapping
    mash_mapping_bin = bin_prefix+".mash_dist.mapping" 
    mash_dist_map_cmd = ["MapPatricGenomeIDs.py","-m",patricMapping,"-f",bin_mash_file,"-o",mash_mapping_bin]
    subprocess.check_call(mash_dist_map_cmd)
    #diamond vfdb
    dmnd_bin_cmd = ["Run_Diamond.py","-t",args.threads,"-f","../"+bin_file,"-d",vfdbDB,"-p",bin_prefix]
    time_dict["Diamond"]["Bin"][bin_prefix]["Start"] = time.time()
    subprocess.check_call(dmnd_bin_cmd)
    time_dict["Diamond"]["Bin"][bin_prefix]["End"] = time.time()
    #TODO: pathogenicity score
    pp_cmd = ["scoreGenomeByVfdbLogodds.py",bin_prefix+".diamond.out"]
    time_dict["PP_Profile"]["Bin"][bin_prefix]["Start"] = time.time()
    pp_bin_stdout = bin_prefix+".pp_profile.out"
    pp_bin_stderr = bin_prefix+".pp_profile.err"
    with open(pp_bin_stdout,"w") as pbo, open(pp_bin_stderr,"w") as pbe:
        subprocess.check_call(pp_cmd,stdout=pbo,stderr=pbe)
    time_dict["PP_Profile"]["Bin"][bin_prefix]["End"] = time.time()

    #TODO: checkm
    #time_dict["Checkm"]["Bin"]["Start"] = time.time()
    #subprocess.check_call(checkm_cmd)
    #time_dict["Checkm"]["Bin"]["End"] = time.time()

    #TODO: anl models

    os.chdir("..")

with open("Print_Time_Dict.txt","w") as ptd:
    ptd.write(json.dumps(time_dict))

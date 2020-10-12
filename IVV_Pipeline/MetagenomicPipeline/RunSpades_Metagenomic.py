import sys, subprocess, argparse, os, shutil, glob

parser = argparse.ArgumentParser()
parser.add_argument("-t","--threads",default="1")
parser.add_argument("-r1","--read1",required=True)
parser.add_argument("-r2","--read2",required=True)
parser.add_argument("-s","--spades_dir",required=True)
parser.add_argument("-p","--prefix",required=True)

args = parser.parse_args()

#add "/" to spades_dir if not present
spades_path = "/scratch/cc8dm/iSentry_scripts/IVV_Pipeline/SPAdes/RunSpades.py"
spades_cmd = [spades_path,"-m","-t",args.threads,"-r",args.r1+","+args.r2,"-o",args.spades_dir]
spades_log_out = args.spades_dir+"log.stdout"
spades_log_err = args.spades_dir+"log.stderr"
with open(spades_log_out,"w") as slo and open(spades_log_err,"w") as sle:
    subprocess.check_call(spades_cmd,stdout=slo,stderr=sle)

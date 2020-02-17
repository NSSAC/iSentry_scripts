import sys
import subprocess
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-t','--threads')
parser.add_argument('-d','--database')
parser.add_argument('-m','--map')
parser.add_argument('-c','--headers')
parser.add_argument('-f','--filename')
parser.add_argument('-s','--sample')

args = parser.parse_args()

outfile = sample+".mash_out.txt"
errfile = sample+".mash_err.txt"

#TODO: Mash theshold list implementation
mash_threshold = 0.0

while not os.path.isfile(outfile) or os.path.getsize(outfile) == 0:
    command = "mash dist -d " + str(mash_threshold) + " -p " + args.threads + " " + args.database + " " + args.filename  

    with open(outfile,"w") as out, open(errfile,"w") as err:
        process = subprocess.Popen(command,shell=True,stdout=out,stderr=err)
        process.wait()

    #TODO: implement iteration over ideal mash thresholds
    #Potentially just a list of values where in index is tracked from the most strict to least strict
    if os.path.getsize(outfile) == 0:
        mash_threshold = mash_threshold * 2

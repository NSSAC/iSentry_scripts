import sys
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-t','--threads')
parser.add_argument('-d','--database')
parser.add_argument('-m','--map')
parser.add_argument('-h','--headers')
parser.add_argument('-f','--filename')

args = parser.parse_args()

mash_threshold = 0.05

command = "mash dist -d " + mash_threshold + " -p " + args.threads + " " + args.database + " " + args.filename  

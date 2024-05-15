#!/usr/bin/env python3

"""
Unpacks consensus sequence files stored in 2bit format on cbsubscb17 to fasta files on danko_0001. Make sure to delete
fasta files when done.
"""

import glob
import os

files = glob.glob("/home/ayh8/CLIPNET/data/gse110638/1000genomes_yrb/consensus/*.2bit")
directory = os.path.split(files[0])[0]

wkdir = "/home/danko_0001/projects/ayh8/1000genomes_yrb/consensus/"
os.chdir(wkdir)

for f in files:
    prefix = os.path.split(f)[1].split(".")[0]
    os.system("cp %s %s" % (f, wkdir))
    os.system("twoBitToFa ./%s.2bit ./%s.fna" % (prefix, prefix))
    os.system("rm ./%s.2bit" % prefix)

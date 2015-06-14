#! /usr/bin/env python
# coding=utf-8

# Wrapper to invoke qucsator with the onelab socket
# set in the modified netlist file.
#
# Alexander Krimm <alex@stud.tu-darmstadt.de>
# TU-Darmstadt, 2015

import subprocess
import sys
import os

name = ""
addr = ""
infile = ""
outfile = ""

# read socket and solver name from arguments
for i, v in enumerate(sys.argv):
    if v == '-onelab':
        name = sys.argv[i + 1]
        addr = sys.argv[i + 2]
    if v == '-i':
        infile = sys.argv[i + 1]
    if v == '-o':
        outfile = sys.argv[i + 1]

# write name and socket in temporary file
tmpfile = os.getcwd() + "/.onelabsocket.tmp"
with open(tmpfile, 'w') as tmpf:
    tmpf.write(name + "\n" + addr)

# call qucsator
subprocess.call(["qucsator", "-b", "-i", infile + ".netlist.tmp",
                 "-o", outfile, "-p", os.getcwd(), "-m", "qucsonelab",
                 "qucsonelab4"])

os.remove(tmpfile)

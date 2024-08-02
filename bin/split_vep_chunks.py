#!/usr/bin/env python3.8

import math
import os
import sys
import subprocess

VEP_FILE = sys.argv[1]
RAMlimit = float(sys.argv[2])

# numHeaderSkip=354#this depends on the type of annotation used
numHeaderSkip = 0
with open(VEP_FILE, "r") as F:
    for line in F:
        if line.startswith("##"):
            numHeaderSkip += 1
        else:
            break

# 06/30/2023: Process the vep annotation result chunk by chunk
command = f'wc -l {VEP_FILE} | cut -d" " -f1'
nlines = (
    int(subprocess.check_output(command, shell=True, text=True).strip())
    - numHeaderSkip
    - 1
)

filesize = os.path.getsize(VEP_FILE)

num_chunk = math.ceil((filesize / 1024 / 1024) / (RAMlimit * 1024 / 23))

chunk_lines = math.ceil(nlines / num_chunk)

lines = []
if num_chunk > 1:
    for i in range(num_chunk - 1):
        lines.append(
            f"{i+1}\t{int(numHeaderSkip+1)}\t{int(numHeaderSkip+2 + i*chunk_lines)}\t{int(numHeaderSkip+1 + (i+1)*chunk_lines)}"
        )


lines.append(
    f"{num_chunk}\t{int(numHeaderSkip+1)}\t{int(numHeaderSkip+2 + (num_chunk-1)*chunk_lines)}\t{int(nlines+numHeaderSkip+1)}"
)

vep_split = open("vep_split.txt", "w")
vep_split.write("\n".join(lines))
vep_split.write("\n")
vep_split.close()

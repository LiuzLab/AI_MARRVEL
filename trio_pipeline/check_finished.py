
from glob import glob
import os
import sys
import subprocess
import pandas as pd
# Check all cases if they are finshed with no error
def check_integraty(out_path):

    results = glob(f"{out_path}/*.trio.csv")

    if len(results) == 0:
        return False

    else:
        nlines = result = subprocess.run(["wc", "-l", results[0]], stdout=subprocess.PIPE, text=True)
        nlines = int(nlines.stdout.split()[0])
        if nlines > 1:
            return True
        else:
            return False



samples = sys.argv[1]
sampleDf = pd.read_csv(samples, sep='\t', header=None)
sampleIDs = sampleDf[0].tolist()

for ID in sampleIDs:
    check = check_integraty(f"output/{ID}/final_matrix")
    if check:
        continue
    else:
        print(f"Run {ID} didn't finish.")






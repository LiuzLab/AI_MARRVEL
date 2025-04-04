#!/usr/bin/env python3.8

import sys

from trio.inh_utils import *

score_file = sys.argv[1]
singleton_rank = sys.argv[2]
inheritance_file = sys.argv[3]
out_path = sys.argv[4]


inh_df = matchID_in_inheritance(score_file, singleton_rank, inheritance_file)
inh_df = matrix_inh(inh_df, singleton_rank)
inh_df.to_csv(out_path, sep="\t")



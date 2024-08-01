#!/usr/env/bin python
import sys

import pandas as pd
import numpy as np
import vamb

"""
Arguments: 1. input file, 2. contigLen used in vamb input, 3. output path
"""


df = pd.read_csv(sys.argv[1], sep = '\t')
df = df[df.contigLen >= int(sys.argv[2])]
bamfiles = [x for x in df.columns if x[-4:]  ==  ".bam"]

abd_matrix = df[bamfiles].to_numpy(dtype = np.float32)
minid = 0.001
refhash = None
abd = vamb.parsebam.Abundance(abd_matrix, bamfiles, minid, refhash)

abd.save(sys.argv[3])


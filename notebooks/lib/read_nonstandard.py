#!/usr/bin/env python
# coding: utf-8

"""Reading especially heavy-weight datasets used in some figures (500 Kb and 1 Mb windows)"""

from . import logger, log_step
from .core import *
from .read_data_basic import *


# Load non-standard epigenetics:

# Load all the stacks 500 Kb into memory
log_step("Collecting epigenetic stacks (500 Kb flank)")
files_stacks_500Kb = glob.glob(
    "/home/agalicina/DANIO/HIC/WD_2021_paper/13_epigenetics/data/stacks_500Kb/*/*.npy"
)
stacks_bins_500Kb = {f.split("/")[-1]: np.load(f) for f in files_stacks_500Kb}
log_step(f"Loaded {len(stacks_bins_500Kb)} stacks (500 Kb flank)")
for f in files_stacks_500Kb:
    k = f.split("/")[-1]
    if k not in stacks_bins_500Kb.keys():
        stacks_bins_500Kb[k] = np.load(f)

del f,k

# Load all the stacks 1 Mb into memory
log_step("Collecting epigenetic stacks (1 Mb flank)")
files_stacks_1Mb = glob.glob(
    "/home/agalicina/DANIO/HIC/WD_2021_paper/13_epigenetics/data/stacks_1Mb/*/*.npy"
)
stacks_bins_1Mb = {f.split("/")[-1]: np.load(f) for f in files_stacks_1Mb}
log_step(f"Loaded {len(stacks_bins_1Mb)} stacks (1 Mb flank)")
for f in files_stacks_1Mb:
    k = f.split("/")[-1]
    if k not in stacks_bins_1Mb.keys():
        stacks_bins_1Mb[k] = np.load(f)

del f,k

# Load non-standard Hi-C stacks:
dct_stacks_hic_500Kb = {}
for source_fountains in ["WT"]:
    log_step(f"Loading 500 Kb snippet stacks for {source_fountains}")
    dct_stacks_hic_500Kb[source_fountains] = fontanka.utils.read_snips(
        f"/home/agalicina/DANIO/HIC/WD_2021_paper/10_fountain_new/data/fountains/{source_fountains}.500Kb.snips.npy"
    )

del source_fountains

#!/usr/bin/env python
# coding: utf-8

from . import logger, log_step
from .core import *
from .read_data_basic import *


# Read all fountains: populate a reusable dict with the list of fountain bases per sample
sources_fountains = ["WT", "Wild-Type_2.75", "TR", "Wild-Type_11", "Wild-Type_25"]

dict_fountains = {}

log_step(f"Loading fountain BEDs for sources: {', '.join(sources_fountains)}")
for source_fountains in sources_fountains:

    # Read fountains:
    df_fount = pd.read_table(
        f"/home/agalicina/DANIO/HIC/WD_2021_paper/10_fountain_new/data/bed/ftop_{source_fountains}_22-03-2022.bed",
        header=None,
        index_col=3,
    )
    df_fount.columns = ["chrom", "start", "end", "SIM"]

    # Store fountains:
    dict_fountains[source_fountains] = df_fount.copy()

    del df_fount

del source_fountains, sources_fountains

# ### Load all snippets
sources_fountains = [
    "WT",  #'WT1', 'WT2',
    "PS",  #'PS1', 'PS2',
    "SN",  #'SN1', 'SN2',
    "NP",  #'NP1', 'NP2',
    "TR",  #'TR1', 'TR2',
    "MZsox19b_5.3",  #'MZsox19b_5.3_1', 'MZsox19b_5.3_2',
    "MZnanog_5.3",  #'MZnanog_5.3_1', 'MZnanog_5.3_2',
    "MZspg_5.3",  #'MZspg_5.3_1', 'MZspg_5.3_2',
    "Wild-Type_2.75",  #'Wild-Type_2.75_1', 'Wild-Type_2.75_2',
    "Wild-Type_5.3",  #'Wild-Type_5.3_1', 'Wild-Type_5.3_2',
    "Wild-Type_11",  #'Wild-Type_11_1', 'Wild-Type_11_2',
    "Wild-Type_25",  #'Wild-Type_25_1', 'Wild-Type_25_2',
    "Wike2021_WT_ZF_Embryos_4hpf_FLAV",
    "Wike2021_WT_ZF_Embryos_4hpf_SGC",
    "Wike2021_WT_ZF_Embryos_4hpf_DMSO",
    "sperm",
]


# Load the dct_stacks_hic of snippets (200 Kb on-diagonal map fragments) for each source
dct_stacks_hic = {}
log_step(f"Loading 200 Kb snippet dct_stacks_hic for sources: {', '.join(sources_fountains)}")
for source_fountains in sources_fountains:
    if source_fountains not in dct_stacks_hic.keys():
        dct_stacks_hic[source_fountains] = fontanka.utils.read_snips(
            f"/home/agalicina/DANIO/HIC/WD_2021_paper/10_fountain_new/data/fountains/{source_fountains}.200Kb.snips.npy"
        ).T # We transpose b/c stackups were created by old cooltools that has the dimensions flipped
del source_fountains


# Generate fountain similarity score per source against the WT template
from fontanka import generate_similarity_score

log_step("Loading WT average fountain template and computing similarity scores")
avfountain = np.load(
    "/home/agalicina/DANIO/HIC/WD_2021_paper/10_fountain_new/data/WT.avfountain.npy"
)
df_fountain_strength = pd.DataFrame()

for source_fountains in sources_fountains:
    if not source_fountains in df_fountain_strength.columns:
        scores = generate_similarity_score(
            dct_stacks_hic[source_fountains], avfountain, measure="corr"
        )
        df_fountain_strength.loc[:, source_fountains] = scores

del source_fountains, scores


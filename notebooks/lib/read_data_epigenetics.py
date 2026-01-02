#!/usr/bin/env python
# coding: utf-8

from . import logger, log_step
from .core import *
from .read_data_basic import *

# ### Stacks of epigenetics
log_step("Collecting epigenetic stacks (default flank)")
files_stacks = glob.glob(
    "/home/agalicina/DANIO/HIC/WD_2021_paper/13_epigenetics/data/stacks/*/*.npy"
)
log_step("Collecting epigenetic background bedgraphs")
files_bgs = glob.glob(
    "/home/agalicina/DANIO/HIC/WD_2021_paper/13_epigenetics/data/bedgraphs/*/*.bg"
)


# Load background bedgraphs for normalizing pileups/profiles
bgs_bins = {
    f.split("/")[-1]: pd.read_table(f, header=None).set_axis(
        ["chrom", "start", "end", "value"], axis=1
    )
    for f in files_bgs
}
log_step(f"Loaded {len(bgs_bins)} background bedgraphs")


# Load all the stacks into memory (default flank)
stacks_bins = {f.split("/")[-1]: np.load(f) for f in files_stacks}
log_step(f"Loaded {len(stacks_bins)} stacks (default flank)")

for f in files_stacks:
    k = f.split("/")[-1]
    if k not in stacks_bins.keys():
        stacks_bins[k] = np.load(f)

del f, k


# ### Dataframe with epigenetics

# Build a genome-wide 10 Kb bin table and add epigenetic tracks
log_step("Building 10 Kb bin annotation table")
df_annotations = binnify(chromsizes, binsize)


# Collect all bedgraph tracks and merge into the annotation table
paths = glob.glob(
    "/home/agalicina/DANIO/HIC/WD_2021_paper/13_epigenetics/data/bedgraphs/*/*.bg"
)


log_step(f"Merging {len(paths)} bedgraph tracks into bin annotations")
for f in tqdm.tqdm(paths):
    title = f.split("/")[-1].rstrip(".bg")
    df_tmp = pd.read_table(f, header=None)
    df_tmp.columns = ["chrom", "start", "end", title]
    df_annotations = pd.merge(
        df_annotations, df_tmp, on=["chrom", "start", "end"], how="left"
    )

del df_tmp, f, title


# ### Load initiation zones


# Load replication initiation zones and derive midpoint coordinates
log_step("Loading replication initiation zones (Supplemental Dataset 4)")
df_IZ = pd.read_excel(
    "/home/agalicina/DANIO/HIC/WD_2021_paper/13_epigenetics/data/tables/Supplemental Dataset 4.xlsx",
    header=1,
)
df_IZ.loc[:, "mid"] = (df_IZ["end"] + df_IZ["start"]) // 2


# Pairwise tables linking origins to fountains (and vice versa) for distance analyses
log_step("Loading origin-to-fountain distance tables")
df_IZs_origin2fount = pd.read_table(
    "/home/agalicina/DANIO/HIC/WD_2021_paper/13_epigenetics/data/IZ_Daria/origins_dome_sorted_by_distance_to_fountains.txt"
)
df_IZs_fount2origin = pd.read_table(
    "/home/agalicina/DANIO/HIC/WD_2021_paper/13_epigenetics/data/IZ_Daria/fountains_sorted_by_distance_to_origins.txt",
    header=None,
)
df_IZs_fount2origin.columns = [
    "chrom",
    "start",
    "end",
    "fountain_bin",
    "number",
    "strand",
    "chrom_IZ",
    "start_IZ",
    "end_IZ",
    "replication_IZ",
    "distance",
]
df_IZs_fount2origin.loc[:, "distance_abs"] = np.abs(
    df_IZs_fount2origin.loc[:, "distance"]
)
df_IZs_origin2fount.loc[:, "mid"] = (
    df_IZs_origin2fount.loc[:, "start"] + df_IZs_origin2fount.loc[:, "end"]
) // 2


base_cols = [
    "chrom",
    "start",
    "end",
    "weight",
    "index_IZ",
    "chr_IZ",
    "start_IZ",
    "end_IZ",
    "mid_IZ",
    "distance to fountain_IZ",
]


# Annotate genome bins with nearest initiation zones and distances
bins_IZ = (
    bioframe.overlap(
        bins,
        df_IZs_origin2fount,
        cols2=["chr", "mid", "mid"],
        return_index=True,
        suffixes=["", "_IZ"],
    )
    .drop_duplicates(subset="index")
    .set_index("index")
    .loc[:, base_cols]
    .sort_index()
)

del base_cols


# ### ATAC-Seq peaks
path = "/home/agalicina/DANIO/HIC/WD_2021_paper/10_fountain_new/data_epigenetics/atac-seq_Daria/WT_ATAC-seq_110_bp_around_summit_dome_oblong_common_peaks.txt"
log_step("Loading ATAC-seq peaks around summits (WT)")
df_atac_peaks = pd.read_table(path)
df_atac_peaks.columns = [
    "chrom",
    "start",
    "end",
    "name",
    "ATAC-seq_peak_score_MACS2",
    "_",
    "_",
    "_",
    "_",
]

del path
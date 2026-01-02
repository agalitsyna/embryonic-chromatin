#!/usr/bin/env python
# coding: utf-8

from . import logger, log_step
from .core import *

# # Load universally used data
binsize = 10_000

log_step("Loading chromsizes for danRer11 (reduced)")
chromsizes = read_chromsizes("~/DANIO/GENOME/danRer11_reduced.chrom.sizes")


log_step("Loading chromosome arm sizes table")
chromarms = pd.read_csv("/home/agalicina/DANIO/GENOME/danRer11.armsizes.txt")
chromarms.loc[:, "name"] = chromarms.apply(
    lambda x: f"{x.chrom}:{x.start}-{x.end}", axis=1
)

# Read the data in mutants:
# This table is produced by ../06_mutants_ann_upd.ipynb; it holds differential fountain calls
log_step("Loading differential fountain calls (mutant vs WT)")
df_fount_differential = pd.read_csv(
    "/home/agalicina/DANIO/HIC/WD_2021_paper/10_fountain_new/differential_fountains.GaussianFit15.02-07-22.csv",
    index_col=0,
)

cols = [x for x in df_fount_differential.columns if "diff" in x]
df_fount_differential.loc[:, cols] = df_fount_differential.loc[:, cols].replace(
    {1: "down", 2: "same", 3: "up"}
)
del cols


# Collect all cooler paths and get their prefixes
log_step("Collecting cooler paths for danRer11 datasets")
PATH = "/home/agalicina/DANIO/HIC/data_danrer11/distiller/results_danrer11/"

# List all coolers:
clr_files = glob.glob(
    f"{PATH}/coolers_library_group/*danrer11-reduced.mapq_30.1000.mcool"
) + glob.glob(f"{PATH}/coolers_library/*danrer11-reduced.mapq_30.1000.mcool")

clr_prefixes = list(
    map(lambda x: x.split(".danrer11-reduced")[0].split("/")[-1], clr_files)
)
dct_files = dict(zip(clr_prefixes, clr_files))



# Load a base cooler and bins table for 10 Kb analyses
binsize = 10000
PATH = "/home/agalicina/DANIO/HIC/data_danrer11/distiller/results_danrer11/coolers_library_group/"
coolpath = f"{PATH}WT.danrer11-reduced.mapq_30.1000.mcool::/resolutions/{binsize}"
log_step(f"Opening base cooler at {coolpath}")
clr = cooler.Cooler(coolpath)
log_step("Extracting bins table from cooler")
bins = clr.bins()[:]

# ### Load clean bins
# Clean bins are those 50 Kb away from any unmapped bin at 10 Kb.
# Load clean bins: regions buffered 50 Kb away from unmappable bins
offset_cov = 50_000
log_step(
    f"Loading clean bins with {offset_cov/1000:.0f} Kb buffer from unmappable regions"
)
df_clean_bins = pd.read_table(
    f"/home/agalicina/DANIO/HIC/WD_2021_paper/05_fountain/data_checkpoint_5Jul2023_fontanka_v0.1/selected.{offset_cov}-safe.danrer11-reduced.tsv"
)


# ### Load CTCF and boundaries
# Read the motifs: JASPAR CTCF hits
log_step("Loading CTCF motif calls from JASPAR track")
df_motifs = pd.read_table(
    "/home/agalicina/CTCF_motifs/JASPAR-UCSC-tracks/tracks/danRer11/MA0139.1_p5e-02_r0.8_bg-a0.317-t0.317-g0.183-c0.183.tsv.gz",
    sep="\t",
    header=None,
)
df_motifs.columns = ["chrom", "start", "end", "name", "score", "p-value", "strand"]


# Read open chromatin regions at dome stage and annotate CTCF motifs inside peaks
prefix = "ATAC-Seq_4.5h_DCD019097DT"
log_step(f"Loading ATAC-seq peaks ({prefix}) and annotating CTCF motifs")
df_peaks = pd.read_table(
    f"/home/agalicina/DANIO/ATACSEQ/DCC/{prefix}.narrowPeak.gz",
    compression="gzip",
    header=None,
)
df_peaks.columns = [
    "chrom",
    "start",
    "end",
    "peak_name",
    "score",
    "_",
    "_1",
    "_2",
    "_3",
    "_4",
]
df_peaks = bioframe.sort_bedframe(df_peaks)

# Annotate ATAC-Seq peaks with a given threshold:
th = 5e-6
df_tmp = df_motifs.loc[df_motifs["p-value"] < th, :]
df_peaks_annotated = bioframe.overlap(
    df_peaks, df_motifs, return_index=False, suffixes=["", "_motif"]
).dropna(subset=["start_motif"])

# Get the list of CTCF motifs at regions of open chromatin:
df_ctcf = bioframe.sort_bedframe(df_peaks_annotated.groupby("peak_name").first())

del th, df_tmp, df_peaks, prefix, df_peaks_annotated


# Read open chromatin regions at 12 hpf and annotate CTCF motifs inside peaks
prefix = "ATAC-Seq_12h_DCD019077DT"
log_step(f"Loading ATAC-seq peaks ({prefix}) and annotating CTCF motifs")
df_peaks_12 = pd.read_table(
    f"/home/agalicina/DANIO/ATACSEQ/DCC/{prefix}.narrowPeak.gz",
    compression="gzip",
    header=None,
)
df_peaks_12.columns = [
    "chrom",
    "start",
    "end",
    "peak_name",
    "score",
    "_",
    "_1",
    "_2",
    "_3",
    "_4",
]
df_peaks_12 = bioframe.sort_bedframe(df_peaks_12)

# Annotate ATAC-Seq peaks with a given threshold:
th = 5e-6
df_tmp = df_motifs.loc[df_motifs["p-value"] < th, :]
df_peaks_12_annotated = bioframe.overlap(
    df_peaks_12, df_motifs, return_index=False, suffixes=["", "_motif"]
).dropna(subset=["start_motif"])

# Get the list of CTCF motifs at regions of open chromatin:
df_ctcf_12 = bioframe.sort_bedframe(df_peaks_12_annotated.groupby("peak_name").first())

del th, df_tmp, df_peaks_12, prefix, df_peaks_12_annotated


log_step("Loading boundary calls (all strengths)")
df_boundaries = pd.read_csv(
    "/home/agalicina/DANIO/HIC/WD_2021_paper/13_epigenetics/data/boundaries/danio-2022/data/boundaries_200Kbp_all/Wild-Type_11.danrer11-reduced.mapq_30.1000.csv"
)
df_boundaries.columns = ["start", "end", "bs_threshold", "window", "chrom", "IS", "BS"]
df_boundaries = (
    bioframe.overlap(df_boundaries, df_clean_bins, how="left", return_index=True)
    .dropna(subset=["index_"])
    .drop(["index", "index_"], axis=1)
)

log_step("Loading strong boundary calls")
df_boundaries_strong = pd.read_csv(
    "/home/agalicina/DANIO/HIC/WD_2021_paper/13_epigenetics/data/boundaries/danio-2022/data/boundaries_200Kbp_strong/Wild-Type_11.danrer11-reduced.mapq_30.1000.csv"
)
df_boundaries_strong.columns = [
    "start",
    "end",
    "bs_threshold",
    "window",
    "chrom",
    "IS",
    "BS",
]
df_boundaries_strong = (
    bioframe.overlap(df_boundaries_strong, df_clean_bins, how="left", return_index=True)
    .dropna(subset=["index_"])
    .drop(["index", "index_"], axis=1)
)


df_boundaries_robust = (
    bioframe.overlap(df_boundaries, df_ctcf_12, how="left", return_index=True)
    .dropna(subset=["index_"])
    .groupby("index")
    .first()
)
df_boundaries_strong_robust = (
    bioframe.overlap(df_boundaries_strong, df_ctcf_12, how="left", return_index=True)
    .dropna(subset=["index_"])
    .groupby("index")
    .first()
)


bins_boundaries_robust = (
    bioframe.overlap(bins, df_boundaries_robust, how="left", return_index=True)
    .sort_values("IS_", ascending=False)
    .dropna(subset=["index_"])
    .index.values
)

print(
    "Number of different types of boundaries:",
    len(df_boundaries),
    len(df_boundaries_strong),
    len(df_boundaries_robust),  # The ones that will be used throughout the notebook
    len(df_boundaries_strong_robust),
)

# Delete those boundaries that we won't use:
del df_boundaries, df_boundaries_strong, df_boundaries_strong_robust


# ### Single dataframe with WT 5.3 fountains:
# Single reference fountains table for WT 5.3
source_fountains = "WT"
log_step("Loading WT fountain reference table")
df_fountains = pd.read_table(
    f"/home/agalicina/DANIO/HIC/WD_2021_paper/10_fountain_new/data/bed/ftop_{source_fountains}_22-03-2022.bed",
    header=None,
    index_col=3,
)
df_fountains.columns = ["chrom", "start", "end", "SIM"]

del source_fountains


N_fountains = len(df_fountains)
print("Total number of fountains: ", N_fountains)


# ### Compartments
# Load; assign A/B compartment labels to fountains based on eigenvector sign
log_step("Loading compartment eigenvector annotations")
df_compartments = pd.read_csv(
    "/home/agalicina/DANIO/HIC/WD_2021_paper/12_supplementaries/data/EV/final_table.withepigenetics-danrer11.25000.RT-phased.csv"
)
df_fountains_withcomp = bioframe.overlap(df_fountains, df_compartments)
df_fountains_withcomp.loc[:, "AB_annotation"] = np.where(
    df_fountains_withcomp["E1_WT_"] > 0, "A", "B"
)
df_fountains_withcomp = (
    df_fountains_withcomp.drop_duplicates(["chrom", "start", "AB_annotation"])
    .groupby(["chrom", "start"])
    .apply(lambda x: " ".join(x["AB_annotation"].sort_values()))
)


# # Some stats on compartments:
# df_fountains_withcomp.reset_index().groupby(0).count() * 100 / N_fountains

from . import logger, log_step

# ### Insulation
# Read insulation table (200 Kb flank) for all bins
log_step("Loading insulation scores (200 Kb flank)")
window_flank = 200_000
df_insulation = pd.read_csv(
    f"/home/agalicina/DANIO/HIC/WD_2021_paper/02_insulation/results/insualtion-table.with-selected.{binsize}-w{window_flank}.csv"
)

del window_flank


# ### Genes and expression
# Load Ensembl gene models and expand attributes
log_step("Loading Ensembl gene models (danRer11)")
df_genes_ens = bioframe.read_table(
    "/home/agalicina/DANIO/GENOME/danRer11.ensGene.gtf", schema="gtf"
)
df_genes_ens = pd.concat(
    [
        df_genes_ens,
        bioframe.sandbox.gtf_io.parse_gtf_attributes(
            df_genes_ens.attributes, kv_sep=" ", item_sep=";", quotechar='"'
        ),
    ],
    axis=1,
)
df_genes_ens.loc[:, "gene_id"] = df_genes_ens.loc[:, "gene_id"].apply(
    lambda x: x.split(".")[0]
)


# Bring in expression (FPKM) matrix and merge with gene models
log_step("Loading expression (FPKM) table and merging into genes")
df_expr = pd.read_csv(
    "/home/agalicina/DANIO/RNASEQ/DATA_EBI/fpkms.tsv", sep="\t", comment="#"
)

df_genes_ens_annotated = pd.merge(
    df_genes_ens, df_expr, left_on="gene_id", right_on="Gene ID"
)
df_genes_ens_annotated = df_genes_ens_annotated.query('feature=="CDS"')


# Derive TSS coordinates (1 bp) from annotated CDS intervals
df_tss_annotated = df_genes_ens_annotated.copy()
df_tss_annotated.loc[:, "start"] = np.where(
    df_tss_annotated.strand == "+",
    df_tss_annotated.loc[:, "start"],
    df_tss_annotated.loc[:, "end"],
)
df_tss_annotated.loc[:, "end"] = df_tss_annotated.loc[:, "start"] + 1


# ### Maternally deposited genes

# Tag maternally deposited transcripts
log_step("Loading maternal transcript annotations")
df_genes_maternal = pd.read_excel(
    "/home/agalicina/DANIO/RNASEQ/DATA_Daria/Maternal_and Zygotic_Transcripts.xlsx"
)
df_genes_maternal.loc[:, "is_maternal"] = True


df_genes_ens_annotated = pd.merge(
    df_genes_ens_annotated,
    df_genes_maternal,
    left_on="transcript_id",
    right_on="ENSEMBL_transcript",
    how="left",
)
df_genes_ens_annotated.loc[:, "is_maternal"] = df_genes_ens_annotated.loc[
    :, "is_maternal"
].replace(np.nan, False)


df_tss_annotated = pd.merge(
    df_tss_annotated,
    df_genes_maternal,
    left_on="transcript_id",
    right_on="ENSEMBL_transcript",
    how="left",
)
df_tss_annotated.loc[:, "is_maternal"] = df_tss_annotated.loc[:, "is_maternal"].replace(
    np.nan, False
)


# ### Read chromatin states

binsize_states = 10000
bins_states = binnify(chromsizes, binsize_states)

df_colors_all = []

for stage in tqdm.tqdm(["Dome", "12hpf", "24hpf"]):
    df_colors = pd.read_table(
        f"/home/agalicina/DANIO/HIC/WD_2021_paper/13_epigenetics/data/states/from_Daria/{stage}_cPADRE.bed",
        header=None,
    )

    df_colors.columns = ["chrom", "start", "end", "state", "_1", "_2", "_3", "_4", "_5"]
    df_colors = df_colors.sort_values(["chrom", "start", "end"]).dropna()
    df_colors.loc[:, "stage"] = stage
    df_colors_all.append(df_colors.copy())

    # Calculate the number of basepairs per 1 Kb bin covered by each state:
    states = np.unique(df_colors.state)
    for s in states:
        cov = bioframe.coverage(bins_states, df_colors.query(f'state=="{s}"'))
        bins_states.loc[:, f"{s}_{stage}"] = cov.coverage

    bins_states.loc[:, f"totalcov_{stage}"] = bioframe.coverage(bins_states, df_colors)[
        "coverage"
    ]

df_colors_all = pd.concat(df_colors_all, axis=0)

del stage, df_colors, states, s, cov
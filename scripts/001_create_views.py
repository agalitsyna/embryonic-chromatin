import bioframe
import pandas as pd

# oryLat:
ignore_chroms = ['chrM']
path_reduced_chromsizes = "./data/genomes/oryLat2/oryLat2.reduced.chromsizes"
path_view_output = "./data/genomes/oryLat2/oryLat.arms.view.tsv"

# Load reduced chromsizes:
chromsizes_reduced = pd.read_table(path_reduced_chromsizes, header=None)
chromsizes_reduced.columns = []

# Fetch chromosome sizes:
chromsizes = bioframe.fetch_chromsizes('oryLat2', as_bed=True).query(f"chrom not in {ignore_chroms}")
chromnames = [x for x in chromsizes_reduced if x in chromsizes.chrom.values]

# Fetch centromeres positions (if known):
cens = bioframe.fetch_centromeres('oryLat2').query(f"chrom not in {ignore_chroms}")

# Make chromosome arms:
arms = bioframe.make_chromarms(chromsizes, cens, cols_chroms=chromsizes.columns)
arms.start = arms.start.astype(int)
arms.end = arms.end.astype(int)

df = pd.DataFrame(clr.chromsizes).reset_index()
df.columns=['chrom', 'end']
df.loc[:, 'start'] = 0
df.loc[:, ['chrom', 'start', 'end']].to_csv("/home/agalicina/DANIO/HIC/data_oryLat2/oryLat2/oryLat.arms.tsv", sep="\t", index=None, header=None)
df = df.set_index('chrom').loc[[x for x in df['chrom'] if 'chr' in x and x!="chrM"], :].reset_index()
df[['chrom', 'start', 'end']].to_csv("/home/agalicina/DANIO/HIC/data_oryLat2/oryLat2/oryLat2.reduced.sizes", sep="\t", index=None, header=False)

arms.loc[:, 'name'] = arms.apply(lambda x: f"{x.chrom}:{x.start}-{x.end}", axis=1)

arms.set_index('chrom').loc[[x for x in df['chrom'] if 'chr' in x and x!="chrM"], :].reset_index().to_csv("/home/agalicina/DANIO/HIC/data_oryLat2/oryLat2/oryLat.arms.view.tsv", sep="\t", index=None, header=None)

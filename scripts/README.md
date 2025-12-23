WARNING: This is legacy code for running 3D genome analysis in embryogenesis of multiple species.

## Input data

Available whole-genome chromatin interactomes in embryogenesis:

| Species                   | Genome    | Data Source           | Data type | 
|---------------------------|-----------|-----------------------|-----------|
| _Danio rerio_             | danRer11  | Wike 2021, Kaaij 2018 | Hi-C      | 
| _Mus musculus_            | mm10      | Du 2017, Ke 2017      | Hi-C      | 
| _Homo sapiens_            | hg38      | Chen 2019             | Hi-C      | 
| _Oryzias latipes_         | oryLat2   | Nakamure 2018         | Hi-C      |
| _Xenopus tropicalis_      | xenTro10  | Liu 2019              | Hi-C      | 
| _Drosophila melanogaster_ | dm6       | Hug 2017              | Hi-C      |  


## Requirements

    conda

At least 8 cores for computations (default number of threads where it possible is set to 8).

## Installation

Create environment `embryonic-chromatin`:
```bash
conda env create -f environment.yml
conda activate embryonic-chromatin
```

Install ipython kernel: 
```bash
python -m ipykernel install --user --name embryonic-chromatin --display-name "Python (embryonic chromatin)"
```

### Uninstall/update notes

Uninstall:
```bash
conda env remove -n embryonic-chromatin
```

Update:
```bash
conda env update --file environment.yml --prune
```

## Run chromatin analysis

### 1. Genome installation

List of genomes: danRer11, mm10, hg38, oryLat2, xenTro10, dm6

#### Light version for pre-downloaded coolers

Generate files: 
- chromosome sizes file
- annotations of genes (download from UCSC, only for _Danio_)

Run: 
```bash
bash scripts/000a_install_genomes_metadata.sh
python scripts/001_create_views.py
```

TODO: add create views script

#### Full version

! Time-consuming

Generate files: 
- fasta with genomic sequences (download from UCSC) 
- bwa index for Hi-C mapping
- chromosome sizes file
- annotations of genes (download from UCSC, only for _Danio_)

Run: 
```bash
bash scripts/000a_install_genomes_metadata.sh
bash scripts/000b_install_genomes.sh
python scripts/001_create_views.py
```

##### Note on _Danio rerio_ genome
Because danRer11 contains chromosomes with heterozygous variant, 
we will remove them and preserve only conventional chromosomes (`danRer11.reduced.chromsizes`).

We also _de novo_ annotated chromosome arms in _Danio_ based on Hi-C maps (`danRer11.armsizes.txt` and `danRer11.centromeres.manual.txt`).

### 2. Hi-C data mapping

#### Light version

Download already mapped coolers from [OSF](https://osf.io/), open storage for scientific datasets. 
You may go to [OSF embryonic-chromatin project](https://osf.io/hg4xs/)
and download data manually, or download automatically by osfclient: 

1. Create OSF account: https://osf.io/

2. Generate OSF token (your profile -> Settings -> Personal access tokens) and copy it.

3. Create OSF config file with credentials: 
```bash
# Content of .osfcli.config

[osf]
username = # your email
project = hg4xs
token = # your token
```

4. List available data and download the coolers:
```bash
osf ls
osf clone ./
mv osfstorage/coolers data/
```

##### Note on the files upload (dev only)

```bash
osf upload -r -U coolers ./
```

#### Full version

! Time and CPU-consuming

Run [distiller-nf version 0.3.3](https://github.com/open2c/distiller-nf/tree/v0.3.3) nextflow pipeline 
for Hi-C data mapping.
Note that distiller has its own dependencies and asks for separate environment.

Distiller supports multiple types of launch with docker, singularity or local run.
Here is an example how to run it locally: 

```bash
git clone https://github.com/open2c/distiller-nf.git
git checkout tags/v0.3.3 -b latest
cd distiller-nf; cp ../configs/project*.yml ./
```

After preparing the environment and setting up the parameters (e.g. `configs/local.local`), for each config file:

```bash
nextflow run distiller.nf -params project_danrer11-reduced.yml
```
Output coolers will be located in `distiller-nf/results*/cooler_library_group/danRer11/*`.
For simplicity, link all the files to the `./data/coolers` folder:

```bash
cp distiller-nf/results_danrer11_reduced/cooler_library_group/* data/coolers/danRer11/
```

For all other config files, repeat the same. 

TODO: add instructions on how to download raw fastqs

### 3. Scaling plots

Scaling plots are based on 1 Kb coolers. The processing is done in the notebook: 

```
./notebooks/003_scalings.ipynb
```

### 4. Calculate insulation and boundaries

Calculate insualtion at 5 Kb resolution and compare insulation at boundaries:

```
./notebooks/004_insulation_and_boundaries.ipynb
```

### 5. Calculate compartments

Calculate compartments at 25 Kb resolution and compartment strength: 

```
./notebooks/005_compartments.ipynb
```


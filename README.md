
## requirements

run with snakemake >7,<8.

Supercell must be installed from github.

```
mamba create -n fca_reanalysis r-base=4 r-remotes \
      bioconductor-scater=1.30.1 bioconductor-scran=1.30.0 \
      bioconductor-scuttle=1.12.0 \
      r-tidyverse=2.0.0
```

```
remotes::install_github("GfellerLab/SuperCell")
```


## Run locally

```
snakemake --use-conda --use-singularity --cores 4 --rerun-incomplete --notemp -kp
```

## Run on Amarel

```
snakemake --profile amarel --use-conda --use-singularity --cores 4 --rerun-incomplete --notemp -kp
```


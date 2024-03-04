
## snakemake

run with snakemake >7,<8.

## Run locally

```
snakemake --use-conda --use-singularity --cores 4 --rerun-incomplete --notemp -kp
```

## Run on Amarel

```
snakemake --profile amarel --use-conda --use-singularity --cores 4 --rerun-incomplete --notemp -kp
```




Use the following command to generate the rulegraph for 1_sbs:

```sh
snakemake --snakefile 1_sbs.smk --rulegraph | dot -Gdpi=300 -Tpng -o 1_sbs_rulegraph.png
```

Use the following command to run the snakemake file:

```sh
snakemake --snakefile 1_sbs.smk
```
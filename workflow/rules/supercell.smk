rule supercell_per_experiment:
    input:
        sce = rules.add_fca_annotation.output.sce
    output:
        sce = "results/supercell/{sample}.usa.fca_annot.supercell.sce.rds",
        res_l = "results/supercell/{sample}.usa.fca_annot.supercell.res_l.rds"
    params:
        gamma = 20,
        k_knn = 5,
    conda:
        "../envs/bioc_singlecell.yaml"
    script:
        "../scripts/supercell_per_experiment.R"
rule supercell_per_experiment:
    input:
        sce = rules.add_fca_annotation.output.sce
    output:
        sce = "results/supercell/{sample}.usa.fca_annot.supercell.sce.rds",
        res_l = "results/supercell/{sample}.usa.fca_annot.supercell.res_l.rds"
    params:
        gamma = 20,
        k_knn = 5,
    #conda:
    #    "../envs/bioc_singlecell.yaml"
    script:
        "../scripts/supercell_all.R"

rule supercell_all:
    input:
        sce = expand("results/import/fca_annotated_r_objs/{sample}.usa.fca_annot.sce.rds", sample = SAMPLES)
    output:
        sce = "results/supercell_all/all.usa.fca_annot.supercell.sce.rds",
        res_l = "results/supercell_all/all.usa.fca_annot.supercell.res_l.rds"
    params:
        gamma = 20,
        k_knn = 5,
    resources:
        time=40,
        mem=128000,
        cpus=1
    #conda:
    #    "../envs/bioc_singlecell.yaml"
    script:
        "../scripts/supercell_all.R"
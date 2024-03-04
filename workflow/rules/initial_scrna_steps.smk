
rule make_r_objs:
    input:
        frydir = rules.alevin_fry_quant.output.dir
    output:
        sce = "results/import/raw_r_objs/{sample}.usa.sce.rds",
    resources:
        time=5,
        mem=20000,
        cpus=2
    conda:
        "../envs/make_r_objs.yaml"
    priority:
        3
    script:
        "../scripts/import_alevin_fry.R"

rule add_fca_annotation:
    input:
        j = lambda wc: config.get("EMBEDDING_JSON").get(pep.get_sample(wc.sample).get("fca_id")),
        a = lambda wc: config.get("ANNOTATION_TSV").get(pep.get_sample(wc.sample).get("fca_id")),
        sce = rules.make_r_objs.output.sce,
    output:
        sce = "results/import/fca_annotated_r_objs/{sample}.usa.fca_annot.sce.rds",
    resources:
        time=5,
        mem=20000,
        cpus=2
    params:
        fca_id = lambda wc: pep.get_sample(wc.sample).get("fca_id")
    conda:
        "../envs/bioc_singlecell.yaml"
    script:
        "../scripts/add_fca_annotation.R"

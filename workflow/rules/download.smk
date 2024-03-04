rule dl_array_express:
    params:
        r1 = lambda wc: pep.get_sample(wc.sample).fq1_uri[int(wc.ss)],
        r2 = lambda wc: pep.get_sample(wc.sample).fq2_uri[int(wc.ss)],
    output:
        r1 = temp("results/fastqs/{sample}_{ss}_r1.fastq.gz"),
        r2 = temp("results/fastqs/{sample}_{ss}_r2.fastq.gz")
    resources:
        time=960,
        mem=20000,
        cpus=2
    threads: 128 # basically ensures only 1 download runs concurrently when running locally  - to ensure arrayexpress doesn't block us
    shell:
        """
        wget -O {output.r1} {params.r1} &&
        wget -O {output.r2} {params.r2}
        """
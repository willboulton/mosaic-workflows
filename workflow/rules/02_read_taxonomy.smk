include: "00_utils.smk"


rule deinterleave:
    input:
        f"{fastq_tmpdir}" + "/{batch}/{sample}-filter-METAG.fastq.gz",
    output:
        fq1=f"{fastq_tmpdir}" + "/{batch}/{sample}-filter-METAG_1.fastq.gz",
        fq2=f"{fastq_tmpdir}" + "/{batch}/{sample}-filter-METAG_2.fastq.gz",
    params:
        batch = set_param('batch')
    conda:
        "bio"
    shell:
        f"cd {fastq_tmpdir}" + """/{params.batch}
        {software.utils}/deinterleave {input}
        """

rule centrifuge_on_contigs:
    input:
        tar=f"{img_folder}" + "/{sample}.tar.gz",
        fna=f"{img_folder}" + "/{sample}.a.fna",
    output:
        f"{taxo_tmpdir}" + "/{batch}/{sample}.a.cf"
    params:
        centrifuge_idx = software.centrifuge_idx
    envmodules:
        "centrifuge/1.0.4"
    shell:
        """
        centrifuge -x {params.centrifuge_idx} 
            -U <(zcat {input.tar} | tar xfO - {input.fna}) 
            -f -p {resources.cpus_per_task} -S {output}
        """

rule centrifuge_on_reads:
    input:
        fq1=f"{fastq_tmpdir}" + "/{batch}/{sample}-filter-METAG_1.fastq.gz",
        fq2=f"{fastq_tmpdir}" + "/{batch}/{sample}-filter-METAG_2.fastq.gz",
    output:
        f"{taxo_tmpdir}" + "/{batch}/{sample}.fastq.cf"
    envmodules:
        "centrifuge/1.0.4"
    shell:
        """
        centrifuge -x {software.centrifuge_idx} 
            -1 {input.fq1}
            -2 {input.fq2}
            -q -p {resources.cpus_per_task} -S {output}
        """

rule mmseqs:
    input:
    output:
        f"{taxo_tmpdir}" + "/{batch}"
    params:
        sample = set_param("sample_id"),
        tmpdir = f"{taxo_tmpdir}" + "/{batch}/{sample}",
        mmseq_db = software.mmseqs_db,
        img_folder = img_folder
    envmodules:
        "MMSeqs2/14"
    shell:
        """
        scripts/run_mmseqs.sh 
            {params.sample}
            {output}
            {params.tmpdir}
            {params.mmseqs_db}
            {resources.cpus_per_task}
            {params.img_folder}
        """
    
rule mmseq_lca:
rule taxometer:
rule aggregate_taxonomy:
rule bwa_mem:
rule split_fastq_and_tidy:


# Output files from this module:
#  - mapping of fastq to fna
#  - taxonomy of reads and contigs
#  - list of unmapped EVU reads
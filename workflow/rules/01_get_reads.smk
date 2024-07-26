rule globus:
    output:
        f"{fastq_tmpdir}" +  "/{sample}-filter-METAG.fastq.gz",
    conda:
        "bio"
    params:
        nersc_fastq_path = set_param('nersc_fastq_path')
    resources:
        globus_workers = 1
    shell:
        """
        TASK_ID=`globus transfer --skip-activation-check \
        {globus.nersc_uuid}:{params.nersc_fastq_path} \
        {globus.ada_uuid}:{params.img_taxon_id}-filter-METAG.fastq.gz`
        T=`echo $TASK_ID | cut -f 3 -d ':'`
        FINISHED=""
        while [ ! "$FINISHED" ]; do
            $FINISHED=`globus task show $T | grep Status | grep SUCCEEDED`
            sleep 60
        done
        mv {home}/{params.img_taxon_id}-filter-METAG.fastq.gz {output}
        """

rule bwa_mem:
    input:
        tar = f"{img_folder}/{{sample}}.tar.gz",
        fq = f"{fastq_tmpdir}/{{batch}}/{{sample}}-filter-METAG.fastq.gz"
    output:
        f"{bam_tmpdir}/{{batch}}/{{sample}}.sorted.bam"
    envmodules:
        "bwa/0.7.17",
        "samtools/1.16.1"
    shell:
        """
        cd {bam_tmpdir}/{wildcards.batch}
        zcat {input.tar} | \\
            tar xfO - {wildcards.sample}/{wildcards.sample}.a.fna > {wildcards.sample}.a.fna
        bwa index -a bwtsw {wildcards.sample}.a.fna
        bwa mem -M -t {threads} \\
            {wildcards.sample}.a.fna {input.fq} | \\
        samtools sort --threads {threads} -o {output} 
        """

rule deinterleave:
    input:
        f"{fastq_tmpdir}/{{batch}}/{{sample}}-filter-METAG.fastq.gz",
    output:
        fq1=f"{fastq_tmpdir}/{{batch}}/{{sample}}-filter-METAG_1.fastq.gz",
        fq2=f"{fastq_tmpdir}/{{batch}}/{{sample}}-filter-METAG_2.fastq.gz",
    conda:
        "bio"
    shell:
        """
        cd {fastq_tmpdir}/{wildcards.batch}
        {software.utils}/deinterleave_parallel {input}
        """

rule get_img_file:

rule get_coverage_file:

rule append_coverage_file:

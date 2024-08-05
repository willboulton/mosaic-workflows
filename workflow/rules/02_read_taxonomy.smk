rule extract_files:
    input:
        f"{img_folder}/{{sample}}.tar.gz"
    output:
        fna = f"{binning_dir}/{{sample}}/{{sample}}.a.fna",
        faa = f"{binning_dir}/{{sample}}/{{sample}}.a.faa",
        phylodist = f"{taxo_tmpdir}/{{sample}}/{{sample}}.a.phylodist.txt",
        map = f"{binning_dir}/{{sample}}/{{sample}}.a.map.txt",
        depth = f"{binning_dir}/{{sample}}/{{sample}}.a.depth.txt"

    shell:
        """
        tar -xzf {input} {wildcards.sample}/{wildcards.sample}.a.fna -O > {output.fna}
        tar -xzf {input} {wildcards.sample}/{wildcards.sample}.a.faa -O > {output.faa}
        tar -xzf {input} {wildcards.sample}/{wildcards.sample}.a.phylodist.txt -O > {output.phylodist}
        tar -xzf {input} {wildcards.sample}/{wildcards.sample}.a.map.txt -O > {output.map}
        tar -xzf {input} {wildcards.sample}/{wildcards.sample}.a.depth.txt -O > {output.depth}
        """


rule centrifuge_on_contigs:
    input:
        fna = f"{binning_dir}/{{sample}}/{{sample}}.a.fna"
    output:
        f"{taxo_tmpdir}/{{sample}}/{{sample}}.a.cf"
    envmodules:
        "centrifuge/1.0.4"
    shell:
        """
        centrifuge -x {software.centrifuge_idx} \\
            -U {input.fna} \\
            -f -p {threads} -S {output}
        """


rule get_unmapped_reads:
    input:
        f"{bam_tmpdir}/{{sample}}/bam/{{sample}}.sorted.bam"
    output:
        f"{fastq_tmpdir}/{{sample}}/{{sample}}-unmapped.fastq.gz"
    envmodules:
        "samtools/1.16.1"
    shell:
        # -f 12 indicates both the read and the mate are unmapped. 
        """
        samtools fastq -@ {threads} -f 12 -0 /dev/null {input} | gzip -c > {output}
        """

rule mmseqs:
    input:
        f"{binning_dir}/{{sample}}/{{sample}}.a.faa"
    output:
        f"{taxo_tmpdir}/{{sample}}/{{sample}}.faa.mms2.tsv"
    params:
        tmpdir = f"{taxo_tmpdir}/{{sample}}/{{sample}}_mmseqs"
    envmodules:
        "MMSeqs2/14"
    shell:
        """
        {software.workflow}/scripts/run_mmseqs.sh \\
            {wildcards.sample} \\
            {taxo_tmpdir}/{wildcards.sample} \\
            {params.tmpdir} \\
            {software.mmseqs_db} \\
            {threads} \\
            {binning_dir}/{wildcards.sample}
        """

rule mmseqs_lca:
    input:
        f"{taxo_tmpdir}/{{sample}}/{{sample}}.faa.mms2.tsv"
    output:
        f"{taxo_tmpdir}/{{sample}}/{{sample}}.faa.mms2lca"
    params:
        taxonomy_files = f"{software.taxdump}"
    conda:
        "bio"
    shell:
        """
        {software.workflow}/scripts/mmseqs_lca.py {input} {params.taxonomy_files} {output}
        """

rule tiara:
    input:
        f"{binning_dir}/{{sample}}/{{sample}}.a.fna"
    output:
        f"{taxo_tmpdir}/{{sample}}/{{sample}}.a.tiara.txt"
    conda:
        "tiara"
    shell:
        """
        tiara -m 2000 -t {threads} -i {input} -o {output}
        """


# Priority is:
#  - first trust taxometer from the phylodist file 
#  - then trust taxometer from mmseqs
#  - then tiara
#  - add any remaining contigs from centrifuge
rule aggregate_taxonomy:
    input:
        phylodist = f"{taxo_tmpdir}/{{sample}}/{{sample}}.a.phylodist.txt",
        mms2 = f"{taxo_tmpdir}/{{sample}}/{{sample}}.faa.mms2lca",
        tiara = f"{taxo_tmpdir}/{{sample}}/{{sample}}.a.tiara.txt",
        cent = f"{taxo_tmpdir}/{{sample}}/{{sample}}.a.cf", 
        map = f"{binning_dir}/{{sample}}/{{sample}}.a.map.txt"
    output:
        f"{taxo_tmpdir}/{{sample}}/{{sample}}.a.euks.txt"
    conda:
        "bio"
    shell: 
        """
        {software.workflow}/scripts/aggregate_taxonomy.py \\
            {input.phylodist} \\
            {input.mms2} \\
            {input.tiara} \\
            {input.cent} \\
            {software.taxdump} \\
            {input.map} \\
            {output}
        """


rule split_fastq:
    input:
        euks = f"{taxo_tmpdir}/{{sample}}/{{sample}}.a.euks.txt",
        bam = f"{bam_tmpdir}/{{sample}}/bam/{{sample}}.sorted.bam"
    output:
        f"{fastq_tmpdir}/{{batch}}/{{sample}}-EVU-mapped.fastq.gz"
    envmodules:
        "samtools/1.16.1"
    shell:
        """
        samtools view -F 12 -@ {threads} --regions-file {input.euks} {input.bam} | \
            samtools fastq -@ {threads} -0 /dev/null | gzip -c > {output}
        """




# Output files from this module:
#  - mapping of fastq to fna
#  - taxonomy of reads and contigs
#  - EUV reads
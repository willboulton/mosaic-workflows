rule extract_files:
    input:
        f"{img_folder}/{{sample}}.tar.gz"
    output:
        fna = f"{binning_dir}/{{sample}}/{{sample}}.a.fna",
        depth = f"{binning_dir}/{{sample}}/{{sample}}.a.depth.txt",
        phylodist = f"{binning_dir}/{{sample}}/{{sample}}.a.phylodist.txt",
        map = f"{binning_dir}/{{sample}}/{{sample}}.a.map.txt",
        faa = f"{binning_dir}/{{sample}}/{{sample}}.a.faa"
    shell:
        """
        tar -xzf {input} > {output.fna}
        tar -xzf {input} > {output.faa}
        tar -xzf {input} > {output.phylodist}
        tar -xzf {input} > {output.depth}
        tar -xzf {input} > {output.map}
        """


rule centrifuge_on_contigs:
    input:
        fna = f"{binning_dir}/{{sample}}/{{sample}}.a.fna"
    output:
        f"{taxo_tmpdir}/{{sample}}.a.cf"
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
        fq1=f"{fastq_tmpdir}/{{sample}}-unmapped_1.fastq.gz",
        fq2=f"{fastq_tmpdir}/{{sample}}-unmapped_2.fastq.gz"
    envmodules:
        "samtools/1.16.1"
    shell:
        # -f 12 indicates both the read and the mate are unmapped. 
        """
        samtools fastq -@ {threads} -f 12 \\
            -1 {output.fq1} \\
            -2 {output.fq2} \\
            -0 /dev/null -s /dev/null -n {input}
        """


rule centrifuge_on_reads:
    input:
        fq1=f"{fastq_tmpdir}/{{sample}}-unmapped_1.fastq.gz",
        fq2=f"{fastq_tmpdir}/{{sample}}-unmapped_2.fastq.gz",
    output:
        all = f"{taxo_tmpdir}/{{sample}}.unmapped.fastq.cf",
        euks = f"{taxo_tmpdir}/{{sample}}/{{sample}}.euks.cf"
    envmodules:
        "centrifuge/1.0.4"
    shell:
        """
        centrifuge -x {software.centrifuge_idx} \\
            -1 {input.fq1} \\
            -2 {input.fq2} \\
            -q -p {threads} -S {output.all}

        grep -P "Eukaryota|Unk|viria" {output.all} | cut -f 1 > {output.euks}
        """

rule mmseqs:
    input:
        f"{binning_dir}/{{sample}}/{{sample}}.a.faa"
    output:
        f"{taxo_tmpdir}/{{sample}}/{{sample}}.faa.mms2"
    params:
        tmpdir = f"{taxo_tmpdir}/{{sample}}/{{sample}}_mmseqs"
    envmodules:
        "MMSeqs2/14"
    shell:
        """
        scripts/run_mmseqs.sh \\
            {wildcards.sample} \\
            {taxo_tmpdir}/{wildcards.batch} \\
            {params.tmpdir} \\
            {software.mmseqs_db} \\
            {threads} \\
            {img_folder}
        """

rule mmseqs_lca:
    input:
        f"{taxo_tmpdir}/{{sample}}/{{sample}}.faa.mms2"
    output:
        f"{taxo_tmpdir}/{{sample}}/{{sample}}.faa.mms2lca"
    params:
        taxonomy_files = f"{software.taxdump}"
    conda:
        "bio"
    shell:
        """
        scripts/mmseqs_lca.py {input} {params.taxonomy_files} {output}
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
        phylodist = f"{binning_dir}/{{sample}}/{{sample}}.a.phylodist.txt",
        mms2 = f"{taxo_tmpdir}/{{sample}}/{{sample}}.faa.mms2lca",
        tiara = f"{taxo_tmpdir}/{{sample}}/{{sample}}.a.tiara.txt",
        cent = f"{taxo_tmpdir}/{{sample}}.unmapped.fastq.cf"
    output:
        f"{taxo_tmpdir}/{{sample}}/{{sample}}.a.euks.txt"
    conda:
        "bio"
    shell: 
        """
        scripts/aggregate_taxonomy.py \\
            {input.phylodist} \\
            {input.mms2} \\
            {input.tiara} \\
            {input.cent} \\
            {output}
        """

rule split_fastq:
    input:
        fq = f"{fastq_tmpdir}/{{sample}}-filter-METAG.fastq.gz",
        bam = f"{bam_tmpdir}/{{sample}}/bam/{{sample}}.sorted.bam",
        unmapped = f"{taxo_tmpdir}/{{sample}}/{{sample}}.euks.cf",
        mapped = f"{taxo_tmpdir}/{{sample}}/{{sample}}.a.euks.txt"
    output:
        mapped = f"{fastq_tmpdir}/{{batch}}/{{sample}}-EVU-mapped.fastq.gz",
        unmapped = f"{fastq_tmpdir}/{{batch}}/{{sample}}-EVU-unmapped.fastq.gz"
    envmodules:
        "bbmap/37.28", 
        "samtools/1.16.1"
    shell:
        """
        samtools fastq -@ {threads} -f 12 {input.bam} | \\
            filterbyname.sh include=true in=stdin.fq.gz names={input.unmapped} out={output.unmapped}
        
        filterbyname.sh include=true in={input.fq} names={input.mapped} out={output.mapped}
        """


# Output files from this module:
#  - mapping of fastq to fna
#  - taxonomy of reads and contigs
#  - EUV reads
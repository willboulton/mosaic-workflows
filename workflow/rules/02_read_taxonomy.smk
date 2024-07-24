#include: "00_utils.smk"


rule deinterleave:
    input:
        f"{fastq_tmpdir}" + "/{batch}/{sample}-filter-METAG.fastq.gz",
    output:
        fq1=f"{fastq_tmpdir}" + "/{batch}/{sample}-filter-METAG_1.fastq.gz",
        fq2=f"{fastq_tmpdir}" + "/{batch}/{sample}-filter-METAG_2.fastq.gz",
    conda:
        "bio"
    shell:
        """
        cd {fastq_tmpdir}/{wildcards.batch}
        {software.utils}/deinterleave_parallel {input}
        """

rule centrifuge_on_contigs:
    input:
        tar = f"{img_folder}" + "/{sample}.tar.gz"
    output:
        f"{taxo_tmpdir}" + "/{batch}/{sample}.a.cf"
    envmodules:
        "centrifuge/1.0.4"
    shell:
        """
        centrifuge -x {software.centrifuge_idx} \\
            -U <(zcat {input.tar} | tar xfO - {wildcards.sample}/{wildcards.sample}.a.fna) \\
            -f -p {resources.cpus_per_task} -S {output}
        """

rule bwa_mem:
    input:
        tar = f"{img_folder}" + "/{sample}.tar.gz",
        fq = f"{fastq_tmpdir}" + "/{batch}/{sample}-filter-METAG.fastq.gz"
    output:
        f"{bam_tmpdir}" + "/{batch}/{sample}.sorted.bam"
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

rule get_unmapped_reads:
    input:
        f"{bam_tmpdir}" + "/{batch}/{sample}.sorted.bam"
    output:
        fq1=f"{fastq_tmpdir}" + "/{batch}/{sample}-unmapped_1.fastq.gz",
        fq2=f"{fastq_tmpdir}" + "/{batch}/{sample}-unmapped_2.fastq.gz"
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
        fq1=f"{fastq_tmpdir}" + "/{batch}/{sample}-unmapped_1.fastq.gz",
        fq2=f"{fastq_tmpdir}" + "/{batch}/{sample}-unmapped_2.fastq.gz",
    output:
        f"{taxo_tmpdir}" + "/{batch}/{sample}.unmapped.fastq.cf"
    envmodules:
        "centrifuge/1.0.4"
    shell:
        """
        centrifuge -x {software.centrifuge_idx} \\
            -1 {input.fq1} \\
            -2 {input.fq2} \\
            -q -p {resources.cpus_per_task} -S {output}
        """

rule mmseqs:
    input:
        f"{img_folder}" + "/{sample}.tar.gz"
    output:
        f"{taxo_tmpdir}" + "/{batch}/{sample}.faa.mms2"
    params:
        tmpdir = f"{taxo_tmpdir}" + "/{batch}/{sample}_mmseqs"
    envmodules:
        "MMSeqs2/14"
    shell:
        """
        scripts/run_mmseqs.sh \\
            {wildcards.sample} \\
            {taxo_tmpdir}/{wildcards.batch} \\
            {params.tmpdir} \\
            {software.mmseqs_db} \\
            {resources.cpus_per_task} \\
            {img_folder}
        """

rule mmseqs_lca:
    input:
        f"{taxo_tmpdir}" + "/{batch}/{sample}.faa.mms2"
    output:
        f"{taxo_tmpdir}" + "/{batch}/{sample}.faa.mms2lca"
    params:
        taxonomy_files = ""
    conda:
        "bio"
    shell:
        """
        scripts/mmseqs_lca.py {input} {params.taxonomy_files} {output}
        """

rule vamb_convert_coverage:
    input:
        f"{bam_tmpdir}" + "/{batch}/{sample}.sorted.bam"
    output:
        jgi_depths="",
        npz=""
    conda:
        "anvio"
    shell:
        """
        jgi_summarize_bam_contig_depths ---
        scripts/jgi_to_vamb_abd.py \\
            {input} \\   # 1. input file
            2000 \\      # 2. contigLen used in vamb input
            {output.npz} # 3. output path
        """

rule taxometer:
    input:
        fna = "",
        bam = ""
    output:
        directory("")
    conda:
        "vamb"
    shell:
        # Perhaps these are better to use - use what's available!
        #  --bamfiles /path/to/bam/*.bam \
        #  --taxonomy another_taxonomy.csv \
        """
        {software}/vamb taxometer \\
            --fasta {input.fna} \\
            --outdir {output} \\
            --bamfiles {input.bam} \\
            --taxonomy {input.phylodist} \\
            --column_contigs contigs \\
            --column_taxonomy taxonomy \\
            --delimiter_taxonomy ,
        """

rule tiara:
    input:
        f"{{binning_dir}}/{sample}.a.fna"
    output:
        f"{{taxo_tmpdir}}/{sample}.tiara.txt"
    conda:
        "tiara"
    threads: 4
    shell:
        """
        tiara -m 2000 -t {threads} -i {input} -o {output}
        """


# Priority is:
#  - first trust taxometer from the phylodist file 
#  - then trust taxometer from mmseqs
#  - add any remaining contigs from centrifuge
rule aggregate_taxonomy:
    input:
        mms2 = "",
        cent = "",
        taxo = ""
    output:
        unmapped = "",
        mapped = ""
    conda:
        "bio"
    script: 
        "aggregate_taxonomy.py {input} {output}"

rule split_fastq:
    input:
        fq = f"{fastq_tmpdir}" + "/{batch}/{sample}-filter-METAG.fastq.gz",
        bam = f"{bam_tmpdir}" + "",
        unmapped = f"{taxo_tmpdir}" + "",
        mapped = f"{taxo_tmpdir}" + ""
    output:
        mapped = f"{fastq_tmpdir}" + "/{batch}/{sample}-EVU-mapped.fastq.gz",
        unmapped = f"{fastq_tmpdir}" + "/{batch}/{sample}-EVU-unmapped.fastq.gz"
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
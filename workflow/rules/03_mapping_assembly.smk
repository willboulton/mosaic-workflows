def collate_fastqs(wc):
    filter_samples = [s for s in samples if s.batch == wc.batch]
    return ({ 
        "paired_mapped": (
            [ f"{fastq_tmpdir}/{s.img_taxon_id}/{s.img_taxon_id}-EVU-paired-mapped.fastq.gz" for s in filter_samples ]         ), 
        "singletons": [ f"{fastq_tmpdir}/{s.img_taxon_id}/{s.img_taxon_id}-EVU-singletons.fastq.gz" for s in filter_samples ]
    })

def collate_fastas(wc):
    filter_samples = [s for s in samples if s.batch == wc.batch]
    return ({
        "fna": [ f"{binning_dir}/{s.img_taxon_id}/{s.img_taxon_id}.a.fna" for s in filter_samples ], 
        "euks": [ f"{taxo_tmpdir}/{s.img_taxon_id}/{s.img_taxon_id}.a.euks.txt" for s in filter_samples ],
        "map": [ f"{binning_dir}/{s.img_taxon_id}/{s.img_taxon_id}.a.map.txt" for s in filter_samples ]
    })


rule concat_fasta:
    input: unpack(collate_fastas)
    output: f"{binning_dir}/{{batch}}/{{batch}}.euks.fna"
    conda: "blobology"
    shell:
        """
        paste \\
            <(echo "{input.euks}" | tr ' ' '\\n') \\
            <(echo "{input.fna}" | tr ' ' '\\n') \\
            <(echo "{input.map}" | tr ' ' '\\n') \\
            <(echo "{input.map}" | tr ' ' '\\n' | sed 's|.a.map.txt||g') | \\
        {software.utils}/parallel -j {threads} -k --colsep '\\t' \\
        "seqkit seq -m 1001 {{2}} | seqkit grep -Inf <(grep -wFf {{1}} {{3}} | cut -f 2) | seqkit replace -p ^ -r {{4/}}_" > {output}
        """


rule concat_fastq:
    input:
        unpack(collate_fastqs)
    output:
        paired_mapped = f"{assembly_dir}/{{batch}}/{{batch}}-EVU.fastq.gz", 
        singletons = f"{assembly_dir}/{{batch}}/{{batch}}-EVU-singletons.fastq.gz"
    shell:
        """
        cat {input.paired_mapped} > {output.paired_mapped}
        cat {input.singletons} > {output.singletons}
        """

rule megahit:
    input:
        interleaved = f"{assembly_dir}/{{batch}}/{{batch}}-EVU.fastq.gz", 
        singletons = f"{assembly_dir}/{{batch}}/{{batch}}-EVU-singletons.fastq.gz"
    output:
        directory(f"{assembly_dir}/{{batch}}/megahit")
    envmodules:
        "megahit/1.2.9"
    shell:
        """
        mkdir -p {assembly_dir}/{wildcards.batch}_tmp
        megahit \\
            --memory 400e9 \\
            --mem-flag 1 \\
            --min-contig-len 500 \\
            --presets meta-large \\
            --tmp-dir {assembly_dir}/{wildcards.batch}_tmp/ \\
            -t {threads} \\
            --12 {input.interleaved} \\
            -r {input.singletons} \\
            --continue -o {output} 
        """

rule strobealign_x_sample:
    input:
        fq = f"{fastq_tmpdir}/{{sample}}/{{sample}}-EVU-paired-mapped.fastq.gz",
        sng = f"{fastq_tmpdir}/{{sample}}/{{sample}}-EVU-singletons.fastq.gz",
        ref = f"{binning_dir}/{{batch}}/{{batch}}.euks.fna"
    output:
        f"{bam_tmpdir}/{{batch}}/bam_multi/{{sample}}_x_{{batch}}_euks_strobalign.sorted.bam"
    conda:
        "anvio"
    envmodules:
        "samtools/1.16.1"
    shell:
        """
        cat {input.fq} {input.sng} > {output}.fastq.gz
        strobealign -t {threads} {input.ref} --interleaved {output}.fastq.gz | \\
            samtools sort -@ {threads} -o {output} && \\
            rm -f {output}.fastq.gz
        """
        # strobealign -t {threads} --aemb {input.ref} --interleaved {input.fq} > {output}

rule strobealign_x_coassembly:
    input:
        fq = f"{fastq_tmpdir}/{{sample}}/{{sample}}-EVU-paired-mapped.fastq.gz",
        sng = f"{fastq_tmpdir}/{{sample}}/{{sample}}-EVU-singletons.fastq.gz",
        ref = f"{assembly_dir}/{{batch}}/megahit/final.contigs.fa"
    output:
        f"{bam_tmpdir}/{{batch}}/bam_coa/{{sample}}_x_{{batch}}_euks_coa.sorted.bam"
    conda:
        "anvio"
    envmodules:
        "samtools/1.16.1"
    shell:
        """
        cat {input.fq} {input.sng} > {output}.fastq.gz
        strobealign -t {threads} {input.ref} --interleaved {output}.fastq.gz | \\
            samtools sort -@ {threads} -o {output} && \\
            rm -f {output}.fastq.gz
        """
        # strobealign -t {threads} --aemb {input.ref} --interleaved {input.fq} > {output}



rule jgi_summarize_bam_x_sample:
    input:
        lambda wc: [ 
            f"{bam_tmpdir}/{s.batch}/bam_multi/{s.img_taxon_id}_x_{s.batch}_euks_strobalign.sorted.bam"
        for s in samples if s.batch == wc.batch ]
    output:
        f"{binning_dir}/{{batch}}/{{batch}}.euks_multi.jgi_cov.tsv"
    conda:
        "anvio"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output} {input}
        """

rule jgi_summarize_bam_coassembly:
    input:
        lambda wc: [ 
            f"{bam_tmpdir}/{s.batch}/bam_coa/{s.img_taxon_id}_x_{s.batch}_euks_coa.sorted.bam"
        for s in samples if s.batch == wc.batch ]
    output:
        f"{binning_dir}/{{batch}}/{{batch}}.euks_coa.jgi_cov.tsv"
    conda:
        "anvio"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output} {input}
        """


rule vamb_coverage_npz_x_sample:
    input:
        f"{binning_dir}/{{batch}}/{{batch}}.euks_multi.jgi_cov.tsv"
    output:
        f"{binning_dir}/{{batch}}/{{batch}}.euks_multi.jgi_cov.npz"
    conda:
        "anvio"
    shell:
        """
        {software.workflow}/scripts/jgi_to_vamb_abd.py {input} 2000 {output}
        """

rule vamb_coverage_npz_coassembly:
    input:
        f"{binning_dir}/{{batch}}/{{batch}}.euks_coa.jgi_cov.tsv"
    output:
        f"{binning_dir}/{{batch}}/{{batch}}.euks_coa.jgi_cov.npz"
    conda:
        "anvio"
    shell:
        """
        {software.workflow}/scripts/jgi_to_vamb_abd.py {input} 2000 {output}
        """


# Output files from this module:
#  - The assembly contigs and scaffolds
#  - The assembly graph
#  - mapping from single sample EVU reads -> assembly
#  - cross sample mappings EVU reads x EVU single assemblies

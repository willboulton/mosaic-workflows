def collate_fastqs(wc):
    filter_samples = [s for s in samples if s.batch == wc.batch]
    return ({ 
        "paired_mapped": (
            [ f"{fastq_tmpdir}/{s.img_taxon_id}/{s.img_taxon_id}-EVU-paired-mapped.fastq.gz" for s in filter_samples ] + 
            [ f"{fastq_tmpdir}/{s.img_taxon_id}/{s.img_taxon_id}-unmapped.fastq.gz" for s in filter_samples ]
        ), 
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

rule kallisto_index:
    input:
        f"{binning_dir}/{{batch}}/{{batch}}.euks.fna"
    output:
        f"{binning_dir}/{{batch}}/{{batch}}.euks.fna.kindex"
    shell:
        """
        {software.kallisto}/kallisto index --distinguish -i {output} -t {threads} {input}
        """

rule deinterleave:
    input:
        f"{fastq_tmpdir}/{{sample}}/{{sample}}-filter-METAG.fastq.gz",
    output:
        fq1=f"{fastq_tmpdir}/{{sample}}/{{sample}}-filter-METAG_1.fastq.gz",
        fq2=f"{fastq_tmpdir}/{{sample}}/{{sample}}-filter-METAG_2.fastq.gz",
    conda:
        "bio"
    shell:
        """
        cd {fastq_tmpdir}/{wildcards.batch}
        {software.utils}/deinterleave_parallel {input}
        """


rule kallisto_x_sample:
    input:
        kindex = f"{binning_dir}/{{batch}}/{{batch}}.euks.fna.kindex", 
        fq1 = f"{fastq_tmpdir}/{{batch}}/{{sample}}-filter-METAG_1.fastq.gz",
        fq2 = f"{fastq_tmpdir}/{{batch}}/{{sample}}-filter-METAG_2.fastq.gz"
    output:
        directory(f"{binning_dir}/{{batch}}/{{sample}}_x_{{batch}}_euks_kallisto")
    shell:
        """
        {software.kallisto}/kallisto quant \\
            -i {input.kindex} \\
            -o {output} \\
            -t {threads} \\
            --plaintext {input.fq1} {input.fq2}
        """

rule strobealign_x_sample:
    input:
        fq = f"{fastq_tmpdir}/{{sample}}-filter-METAG.fastq.gz",
        ref = f"{binning_dir}/{{batch}}/{{batch}}.euks.fna"
    output:
        f"{binning_dir}/{{batch}}/{{sample}}_x_{{batch}}_euks_strobalign.tsv"
    conda:
        "anvio"
    shell:
        """
        strobealign -t {threads} --aemb {input.ref} --interleaved {input.fq} > {output}
        """
        

# Output files from this module:
#  - The assembly contigs and scaffolds
#  - The assembly graph
#  - mapping from single sample EVU reads -> assembly
#  - cross sample mappings EVU reads x EVU single assemblies
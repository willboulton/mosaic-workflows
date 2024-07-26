def collate_fastqs(wc):
    filter_samples = [s for s in samples if s.batch in run_batch]
    return (
        [   f"{fastq_tmpdir}/{s.batch}/{s.img_taxon_id}-EVU-mapped.fastq.gz" for s in filter_samples ] + 
        [ f"{fastq_tmpdir}/{s.batch}/{s.img_taxon_id}-EVU-unmapped.fastq.gz" for s in filter_samples ]
    )

rule concat_fastq:
    input:
        collate_fastqs
    output:
        f"{assembly_dir}/{{batch}}/{{batch}}-EVU.fastq.gz"
    shell:
        """
        cat {input} > {output}
        """

rule megahit:
    input:
        f"{assembly_dir}/{{batch}}/{{batch}}-EVU.fastq.gz"
    output:
        directory(f"{assembly_dir}/{{batch}}/megahit")
    envmodules:
        "megahit/1.2.9"
    shell:
        """
        megahit \\
            --min-contig-length 500 \\
            --presets meta-large \\
            --tmp-dir {assembly_dir}/{wildcards.batch}_tmp \\
            -t {threads} \\
            --12 {input} \\
            -o {output}
        """

rule kallisto_coassembly:

rule kallisto_x_sample:
    
rule assembly_checkpoint_tidyup:

# Output files from this module:
#  - The assembly contigs and scaffolds
#  - The assembly graph
#  - mapping from single sample EVU reads -> assembly
#  - cross sample mappings EVU reads x EVU single assemblies
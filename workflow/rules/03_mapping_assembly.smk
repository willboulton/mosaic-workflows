def collate_fastqs(wc):
    filter_samples = [s for s in samples if s.batch in run_batch]
    return ({ 
        "paired_mapped": (
            [ f"{fastq_tmpdir}/{s.img_taxon_id}/{s.img_taxon_id}-EVU-paired-mapped.fastq.gz" for s in filter_samples ] + 
            [ f"{fastq_tmpdir}/{s.img_taxon_id}/{s.img_taxon_id}-unmapped.fastq.gz" for s in filter_samples ]
        ), 
        "singletons": [ f"{fastq_tmpdir}/{s.img_taxon_id}/{s.img_taxon_id}-EVU-singletons.fastq.gz" for s in filter_samples ]
    })

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
            --min-contig-len 500 \\
            --presets meta-large \\
            --tmp-dir {assembly_dir}/{wildcards.batch}_tmp \\
            -t {threads} \\
            --12 {input.interleaved} \\
            -r {input.singletons} \\
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
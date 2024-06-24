include: "00_utils.smk"

rule concat_fastq:
    input:
        directory(f"{fastq_tmpdir}" + "/{batch}")
    output:
        f"{assembly_dir}" + "/{batch}/{batch}-EVU.fastq.gz"
    shell:
        """
        cat {input}/*-EVU-*.fastq.gz > {output}
        """

rule megahit:
    input:
        f"{assembly_dir}" + "/{batch}/{batch}-EVU.fastq.gz"
    output:
        directory(f"{assembly_dir}" + "/{batch}/megahit")
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

rule bwa_coassembly:

rule bwa_x_sample:
    
rule assembly_checkpoint_tidyup:

# Output files from this module:
#  - The assembly contigs and scaffolds
#  - The assembly graph
#  - mapping from single sample EVU reads -> assembly
#  - cross sample mappings EVU reads x EVU single assemblies
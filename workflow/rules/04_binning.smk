include: "00_utils.smk"

rule vamb_single_assembly:
rule vamb_coassembly:
rule vamb_multiassembly:
rule metabat_single_assembly:
    input:
        img_file = f"{img_folder}" + "/{sample}.tar.gz",
        fna = "{sample}/{sample}.a.fna",
        depth = "{sample}/{sample}.a.depth.txt"
    output:
        directory(f"{binning_dir}/{sample}/metabat")
    conda:
        "anvio"
    shell:
        """
        cd {output}
        tar xzf {input.img_file} {input.fna}
        tar xzf {input.img_file} {input.depth}
        sort -k 2 -t '_' -n {input.depth} > {input.depth}.sorted
        metabat2 \\
            -i {input.fna} \\
            --cvExt -a {input.depth}.sorted \\
            -o mbat. -t {threads}
        rm -rf {sample}
        """

rule metabat_coassembly:
rule vacuum_single_assembly:
rule vacuum_coassembly:
rule vacuum_multiassembly:
rule semibin_single_assembly:
rule semibin_coassembly:
rule semibin_multiassembly:
rule metawrap:
rule euk_bin_size_filtering:
rule virus_sorting:
rule virsorter2:
